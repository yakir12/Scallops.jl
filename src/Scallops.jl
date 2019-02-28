module Scallops

using RayTraceEllipsoids, HCubature, Statistics, Optim, StaticArrays, DataFrames, CSV, Interpolations
export createobj, raytrace_sensitivity, main_figure
RETINAS = (:dis_ret, :pro_ret)
mutable struct Membrane
    passthrough::Float64
    positive::Bool
    rxy::Float64
    rz::Float64
    cz::Float64
    θ0::Float64
    Membrane(passthrough::Float64, positive::Bool, rxy::Real, rz::Real) = new(passthrough, positive, rxy, rz, 0.0, 0.0)
end
Membrane(p::Float64, d::Bool, r::Real) = Membrane(p, d, r, r)
Membrane() = Membrane(1.0, false, 0)
mutable struct Tissue
    ri::Float64
    thick::Float64
    absorption_coefficient::Float64
end
Tissue() = Tissue(0.0, 0.0, 0.0)
struct Layer
    dis_membrane::Membrane
    tissue::Tissue
    pro_membrane::Membrane
end
Layer() = Layer(Membrane(), Tissue(), Membrane())
function getlayers(absorption_coefficient, water_ac, reflectance)
    water_ri = 1.334
    sky              = Membrane(1., true, Inf, Inf)
    water            = Tissue(water_ri, 0, water_ac)
    dis_cornea       = Membrane(1., true, 244, 227)
    cornea           = Tissue(1.37, 23, water_ac)
    dis_lens         = Membrane(1., true, 225, 213)
    lens             = Tissue(1.42, 215, water_ac)
    dis_dis_ret      = Membrane(1., false, 337 - 6)
    dis_ret          = Tissue(1.35, 12, absorption_coefficient[:dis_ret])
    pro_dis_ret      = Membrane(1., false, 337 + 6)
    dis_ret2prox_ret = Tissue(1.34, 102-12-30, water_ac)
    dis_pro_ret      = Membrane(1., false, 337 - 15)
    pro_ret          = Tissue(1.35, 30, absorption_coefficient[:pro_ret])
    pro_pro_ret      = Membrane(1., false, 337 + 15)
    pro_ret2mirror   = Tissue(1.34, 136-102, water_ac)
    mirror           = Membrane(reflectance, false, 417)
    behind_mirror    = Tissue(0, Inf, Inf)
    earth            = Membrane(1., false, Inf)
    (water = Layer(sky, water, dis_cornea),
     cornea = Layer(dis_cornea, cornea, dis_lens), 
     lens = Layer(dis_lens, lens, dis_dis_ret), 
     dis_ret = Layer(dis_dis_ret, dis_ret, pro_dis_ret), 
     dis_ret2prox_ret = Layer(pro_dis_ret, dis_ret2prox_ret, dis_pro_ret), 
     pro_ret = Layer(dis_pro_ret, pro_ret, pro_pro_ret), 
     pro_ret2mirror = Layer(pro_pro_ret, pro_ret2mirror, mirror), 
     mirror = Layer(mirror, behind_mirror, earth))
end
function norefraction!(layers)
    ks = collect(keys(layers))
    todo = setdiff(ks, [:water, :mirror])
    ri = layers.water.tissue.ri
    for i in todo
        layers[i].tissue.ri = ri
    end
    layers
end
function morphlayer!(l, morphz, morphxy)
    l.dis_membrane.rxy *= morphxy
    l.dis_membrane.rz *= morphz
    l.tissue.thick *= morphz
end
function morph2layers(morphz, absorption_coefficient, water_ac, reflectance)
    morphxy = sqrt(1/morphz)
    layers = getlayers(absorption_coefficient, water_ac, reflectance)
    for l in layers
        morphlayer!(l, morphz, morphxy)
    end
    z = sum(l.tissue.thick for l in layers if !isinf(l.tissue.thick))/2
    for l in layers
        l.dis_membrane.cz =  z - signit(l.dis_membrane.positive)*l.dis_membrane.rz
        z -= l.tissue.thick
    end
    layers
end
signit(positive::Bool) = -(-1)^positive
x2h(x, a, b, positive) = signit(positive)*b/a*sqrt(a^2 - x^2)
getzs(l::Layer, photoreceptor_radius) = map([:pro_membrane, :dis_membrane]) do membrane
    m = getfield(l, membrane)
    h = x2h(photoreceptor_radius, m.rxy, m.rz, m.positive)
    m.cz + h
end
function _getmedium(l, photoreceptor_radius)
    z1, z2 = getzs(l, photoreceptor_radius)
    Retina(Cylinder(photoreceptor_radius, z1, z2), l.tissue.absorption_coefficient)
end
getmedium(::Val{:dis_ret}, l, photoreceptor_radius) = _getmedium(l, photoreceptor_radius)
getmedium(::Val{:pro_ret}, l, photoreceptor_radius) = _getmedium(l, photoreceptor_radius)
getmedium(_, l, _) = NonRetina(l.tissue.absorption_coefficient) 
function getopticunits(layers, photoreceptor_radius, ignore_incoming)
    nlayers = length(layers)
    mediums1 = [getmedium(Val(k), v, photoreceptor_radius) for (k,v) in pairs(layers)]
    iₙ = findfirst(isequal(:dis_ret), keys(layers))
    medium_inds = [1:nlayers - 1; nlayers - 1:-1:iₙ]
    ns_ind1 = [1:nlayers - 1; nlayers - 1:-1:iₙ]
    ns_ind2 = [2:nlayers; nlayers - 2:-1:iₙ - 1]
    pointin_inds = [2:nlayers; nlayers-1:-1:iₙ]
    funs = [fill(!, nlayers-1); fill(identity, nlayers - iₙ)]
    surface_inds = [2:nlayers; nlayers - 1:-1:iₙ]
    mediums = mediums1[medium_inds]
    if ignore_incoming
        for retina in RETINAS
            i = findfirst(isequal(retina), keys(layers))
            mediums[i] = NonRetina(layers[1].tissue.absorption_coefficient)
        end
    end
    nss = [layers[i1].tissue.ri/layers[i2].tissue.ri for (i1, i2) in zip(ns_ind1, ns_ind2)]
    pointins = [fun(layers[i].dis_membrane.positive) for (i, fun) in zip(pointin_inds, funs)]
    ellipsoids = [cspheroid(layers[i].dis_membrane.cz, layers[i].dis_membrane.rxy, layers[i].dis_membrane.rz, layers[i].dis_membrane.positive) for i in surface_inds]
    ks = keys(layers)
    names = [i < nlayers ? Symbol(ks[j], :_in) : ks[j] for (i,j) in enumerate(medium_inds)]
    passes = [layers[i].dis_membrane.passthrough for i in surface_inds]
    NamedTuple{Tuple(names)}(OpticUnit.(ellipsoids, pointins, nss, mediums, passes))
end
angulardiff(r::TraceRay, θ) = abs(acos(r.dir[3]) - abs(θ))
angulardiff(r::DeadRay, θ) = abs(θ)
function findnodalz(lb, ub, ous, rc)
    θs = range(-0.02, stop=0.1, length = 4)
    function fun(nodalz, θ)
        l = Light(10^4, rc.rxy, rc, θ, nodalz)
        r = l(TraceRay, 0.5)
        r = raytrace!(ous, r)
        angulardiff(r, θ)
    end
    nodals = map(θs) do θ
        res = optimize(nodalz -> fun(nodalz, θ), lb, ub)
        @assert Optim.minimum(res) < 1e-4 "failed to find optimum for the nodal point"
        # Optim.minimum(res) > 1e-4 && println("failed nodalz")
        Optim.minimizer(res)
    end
    itp = interpolate(nodals, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, θs)
    sitp(0)
end
function createobj(morphz, photoreceptor_radius, absorption_coefficient; ignore_incoming=false, water_ac = 0.03e-6, reflectance = 0.9)
    layers = morph2layers(morphz, absorption_coefficient, water_ac, reflectance)
    # norefraction!(layers)
    mem = layers[1].pro_membrane
    rc = (rxy = mem.rxy, rz = mem.rz, cz = mem.cz)
    ous = getopticunits(layers, photoreceptor_radius, ignore_incoming)
    nodalz = findnodalz(layers.mirror.dis_membrane.cz - layers.mirror.dis_membrane.rz, layers.mirror.dis_membrane.cz + layers.mirror.dis_membrane.rz, ous, rc)
    fls = [abs(mean((layers[retina].dis_membrane.cz - layers[retina].dis_membrane.rz, layers[retina].pro_membrane.cz - layers[retina].pro_membrane.rz)) - nodalz) for retina in RETINAS]
    focallengths = NamedTuple{RETINAS}(fls)
    (layers = layers, rc = rc, ous = ous, nodalz = nodalz, focallengths = focallengths)
end
function raytrace_fraction(distance, aperture, rc, ous, nodalz)
    l = Light(distance, aperture, rc, 0, nodalz)
    function fun(b)
        for retina in RETINAS
            ous[retina].medium.signal.photoreceptor = 0.0
            ous[retina].medium.signal.retina = 0.0
        end
        r = l(Ray, b[1], b[2])
        raytrace!(ous, r)
        SVector(ous[RETINAS[1]].medium.signal.photoreceptor, ous[RETINAS[1]].medium.signal.retina, ous[RETINAS[2]].medium.signal.photoreceptor, ous[RETINAS[2]].medium.signal.retina)
    end
    s, _ = hcubature(fun, SVector(0., 0.), SVector(1., 0.5), initdiv=20)#, maxevals=10^4)
    NamedTuple{RETINAS}(((photoreceptor = 2s[1], retina = 2s[2]), (photoreceptor = 2s[3], retina = 2s[4])))
end
function height2sigma(signal)
    h = signal.photoreceptor/signal.retina
    1/(sqrt(2π)*h)
end
sigma2fwhm(σ) = 2sqrt(2log(2))*σ
spatial2angular(fwhm, focallength) = 2atan(fwhm/2/focallength)
sigma2angular(σ, focallength) = spatial2angular(sigma2fwhm(σ), focallength)
height2angular(signal, focallength) = sigma2angular(height2sigma(signal), focallength)
function fwhm_height(distance, aperture, rc, ous, nodalz, focallengths)
    ss = raytrace_fraction(distance, aperture, rc, ous, nodalz)
    NamedTuple{RETINAS}(([height2angular(ss[retina], focallengths[retina]) for retina in RETINAS]))
end
function optimize_morth4fwhm(retina, distance, aperture, photoreceptor_radius, absorption_coefficient)
    function fun(morphz)
        layers, rc, ous, nodalz, focallengths = createobj(morphz, photoreceptor_radius, absorption_coefficient)
        fwhm = fwhm_height(distance, aperture, rc, ous, nodalz, focallengths)
        fwhm[retina]
    end
    res = optimize(fun, 1, 1.5)#, Brent(), abs_tol = 1e-2)
    Optim.minimizer(res)#, Optim.minimum(res)
end
function optimize_morth(retina, distance, aperture, photoreceptor_radius, absorption_coefficient)
    function fun(morphz)
        layers, rc, ous, nodalz, focallengths = createobj(morphz, photoreceptor_radius, absorption_coefficient)
        Δθ = getΔΘ(focallengths[retina], photoreceptor_radius)
        cosa = cos(Δθ/2)
        s = raytrace_fraction(distance, aperture, rc, ous, nodalz, cosa, retina)
        1/s
    end
    res = optimize(fun, 1, 1.5, Brent(), abs_tol = 1e-2)
    Optim.minimizer(res)#, Optim.minimum(res)
end
getΔΘ(focallength, photoreceptor_radius) = 2atan(photoreceptor_radius/focallength)
modified_photons(D, Δθ) = (π/4)^2*D^2*Δθ^2
modified_land(D, Δθ, k, l) = modified_photons(D, Δθ)*RayTraceEllipsoids.absorption(l, k)
function main_figure(; min_aperture = 200, max_aperture = 400, n_apertures = 8, distance = 10^5, photoreceptor_radius = 5, absorption_coefficient = (dis_ret = 0.00667, pro_ret = 0.00667), filename = joinpath(@__DIR__(), "data"))
    apertures = range(min_aperture, stop = max_aperture, length = n_apertures)
    morphz = 1.1
    layers, rc, ous, nodalz, focallengths = createobj(morphz, photoreceptor_radius, absorption_coefficient)
    df = DataFrame(retina = Symbol[], morpah_factor = Float64[], aperture = Float64[], focallength = Float64[], Δθ = Float64[], dis_ret_ac = Float64[], pro_ret_ac = Float64[], outer_segment_length = Float64[], theoretical_sensitivity = Float64[], raytrace_sensitivity = Float64[], fwhm = Float64[])
    allowmissing!(df)
    for (i, aperture) in enumerate(apertures)
        fwhm = fwhm_height(distance, aperture, rc, ous, nodalz, focallengths)
        for (j, retina) in enumerate(RETINAS)
            Δθ = getΔΘ(focallengths[retina], photoreceptor_radius)
            l = layers[retina].tissue.thick
            theo_sen = modified_land(aperture, δθ, absorption_coefficient[retina], l)
            raytrace_sen = raytrace_sensitivity(distance, aperture, rc, nodalz, ous, retina)
            push!(df, (retina, morphz, aperture, focallengths[retina], Δθ, absorption_coefficient[:dis_ret], absorption_coefficient[:pro_ret], l, theo_sen, raytrace_sen, rad2deg(fwhm[retina])))
        end
    end
    CSV.write("$(filename)dis.csv", df[df[:retina] .== :dis_ret,2:end])
    CSV.write("$(filename)pro.csv", df[df[:retina] .== :pro_ret,2:end])
    CSV.write("$filename.csv", df)
end
function burn!(x)
    for i in length(x)-1:-1:1
        if x[i] < x[i+1]
            x[i] = x[i+1]
        end
    end
end
function getΔΘ(distance, aperture, rc, nodalz, ous, retina)
    function θ2signal(θ)
        l = Light(distance, aperture, rc, θ, nodalz)
        ous[retina].medium.signal.photoreceptor = 0.0
        count = raytrace!(ous, l)
        ous[retina].medium.signal.photoreceptor/count
    end
    maxθ = retina == :dis_ret ? 0.2 : 0.6
    x = range(0, maxθ, length=25)
    y = θ2signal.(x)
    burn!(y)
    itp = interpolate(y, BSpline(Quadratic(Reflect(OnCell()))))
    sitp = scale(itp, x)
    x2 = range(0, maxθ, length = 1000)
    hm = sitp(0)/2
    for x in x2
        if sitp(x) < hm
            return 2x
        end
    end
end
function raytrace_fraction(distance, aperture, rc, ous, nodalz, cosa, retina)
    count = 0
    ous[retina].medium.signal.photoreceptor = 0.0
    n = 1000
    for b in range(0, 1, length=n)
        θ = π - acos(b*(cosa - 1) - cosa)
        l = Light(distance, aperture, rc, θ, nodalz)
        count += raytrace!(ous, l, n = n)
    end
    ous[retina].medium.signal.photoreceptor/count
end
function unmodified_photons(distance, cosa, aperture, rc, nodalz)
    Sₑ = 2π*distance^2*(1 - cosa)
    Sₐ = π*(aperture/2)^2
    h = rc.rz*sqrt(1 - (aperture/2/rc.rxy)^2) + rc.cz
    D = distance + nodalz - h
    Sₑ*Sₐ/D^2
end
function raytrace_sensitivity(distance, aperture, rc, nodalz, ous, retina)
    Δθ = getΔΘ(distance, aperture, rc, nodalz, ous, retina)
    cosa = cos(Δθ/2)
    s = raytrace_fraction(distance, aperture, rc, ous, nodalz, cosa, retina)
    Fₐ = unmodified_photons(distance, cosa, aperture, rc, nodalz)
    s*Fₐ
end

end



