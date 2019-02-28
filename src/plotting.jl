##plottting
using Statistics, Scallops, RayTraceEllipsoids
using Makie, Colors, GeometryTypes
using AbstractPlotting: textslider
x2h(x, a, b, positive) = Scallops.signit(positive)*b/a*sqrt(a^2 - x^2)
x2y(x, a, b, c, positive) = c + x2h(x, a, b, positive)
function Δx2y(x, m1, m2) 
    abs(x2y(x, m1.rxy, m1.rz, m1.cz, m1.positive) - x2y(x, m2.rxy, m2.rz, m2.cz, m2.positive))
end
m1ltm2(x, m1, m2) = x2y(x, m1.rxy, m1.rz, m1.cz, m1.positive) < x2y(x, m2.rxy, m2.rz, m2.cz, m2.positive)
function getintersection(m1, m2)
    r = min(m1.rxy, m2.rxy)
    xs = r/2:r
    i = findfirst(x -> m1ltm2(x, m1, m2), xs)
    xs[i - 1]
end
y2h(y, c) = y - c
y2x(y, a, b, c, positive) = Scallops.signit(positive)*a*sqrt(1 - ((y - c)/b)^2)
function populateθ0!(layers)
    lensl = layers[:lens]
    x = getintersection(lensl.dis_membrane, lensl.pro_membrane)
    m = lensl.dis_membrane
    y = x2y(x, m.rxy, m.rz, m.cz, m.positive)
    for l in layers
        m = l.dis_membrane
        m.θ0 = atan(y2h(y, m.cz), y2x(y, m.rxy, m.rz, m.cz, m.positive))
    end
    layers
end
function θ2coordinate(θ, a, b, c)
    sinθ, cosθ = sincos(θ)
    r = a*b/sqrt((b*cosθ)^2 + (a*sinθ)^2)
    r*cosθ, c + r*sinθ
end
function ellipse_coordinates!(points, θ0, a, b, c, n, keepcorners)
    if θ0 < zero(θ0)
        if keepcorners
            θs = range(π - θ0, stop=2π + θ0, length = n)
        else
            θs = range(π - θ0, stop=2π + θ0, length = n + 1)[2:end]
        end
    else
        if keepcorners
            θs = range(θ0, stop=π - θ0, length = n)
        else
            θs = range(θ0, stop=π - θ0, length = n + 1)[1:end-1]
        end
    end
    for i in 1:n
        points[i] = θ2coordinate(θs[i], a, b, c)
    end
end
function getpoly(layer, n)
    points = Array{Point2f0}(undef, 2n)
    m1 = layer.dis_membrane
    m2 = layer.pro_membrane
    keepcorners = sign(m1.θ0) == sign(m2.θ0)
    ellipse_coordinates!(view(points, 1:n), m1.θ0, m1.rxy, m1.rz, m1.cz, n, keepcorners)
    ellipse_coordinates!(view(points, 2n:-1:n+1), m2.θ0, m2.rxy, m2.rz, m2.cz, n, keepcorners)
    points
end
function morph2layers2(morphz)
    layers = morph2layers(morphz)
    populateθ0!(layers)
end
function getvertices(layers, k, n)
    l = layers[k]
    getpoly(l, n)
end
layers2color(layers, k) = Gray(1 - (layers[k].tissue.ri - 1.32)/0.18)
layers2color(layers, k, m, M) = Gray(1 - (layers[k].tissue.ri - m)/(M - m))
struct ColoredPoint
    coord::Point2f0
    color::Gray{Float64}
end
ColoredPoint(r::DeadRay) = ColoredPoint(Point2f0(NaN, NaN), 1.)
ColoredPoint(r::AbstractRay) = ColoredPoint(Point2f0(r.orig[1], r.orig[3]), 1 - r.int)
function fillpoints!(points, l, b, ous)
    i = 1
    r = l(EmptyRay, b)
    points[i] = ColoredPoint(r)
    for ou in ous
        r = trace!(ou, r)
        i += 1
        points[i] = ColoredPoint(r)
    end
end
function getrays(ous, l; nrays = 21, ϵ = 1e-4)
    npoints = length(ous) + 1
    rays = [[ColoredPoint((NaN, NaN), Gray(0)) for j in 1:npoints] for i in 1:nrays]
    for (ray, b) in zip(rays, range(ϵ, stop=1 - ϵ, length=nrays))
        fillpoints!(ray, l, b, ous)
    end
    rays
end
mm2μm(x) = x*1000
function getranges(min_aperture, max_aperture, aperture_n, min_morph, max_morph, morph_n, min_distance, max_distance, distance_n, min_θ, max_θ, θ_n)
    apertures = range(min_aperture, stop = max_aperture, length = aperture_n)
    morphzs = collect(range(min_morph, stop = max_morph, length = morph_n - 1))
    push!(morphzs, 1)
    sort!(morphzs)
    distances = exp10.(range(log10(min_distance), stop = log10(max_distance), length = distance_n)) #in mm
    θs = range(min_θ, stop = max_θ, length = θ_n)
    (aperture = apertures, morphz = morphzs, distance = distances, θ = θs)
end
function getobservables(apertures, morphzs, distances, θs)
    aperture_slider, aperture = textslider(apertures, "Aperture (μm)", start=mean(apertures))
    morphz_slider, morphz = textslider(morphzs, "Morph", start=1)
    distance_slider, distance = textslider(distances, "Distance (mm)")
    θ_slider, θ = textslider(θs, "Angle (°)", start=0)
    (aperture = aperture_slider, morphz = morphz_slider, distance = distance_slider, θ = θ_slider), (aperture = aperture, morphz = morphz, distance = distance, θ = θ)
end
function ray_trace_applet(; min_aperture = 200, max_aperture = 400, aperture_n = max_aperture - min_aperture + 1, min_morph = 0.8, max_morph = 1.48, morph_n = 100, min_distance = 0.5, max_distance = 1000, distance_n = 100, min_θ = -38, max_θ = 38, θ_n = max_θ - min_θ + 1, limits = FRect(-400,-400,800,800), nvertices = 11, photoreceptor_radius = 5, absorption_coefficient = (dis_ret = 0.00667, pro_ret = 0.00667))
    ranges = getranges(min_aperture, max_aperture, aperture_n, min_morph, max_morph, morph_n, min_distance, max_distance, distance_n, min_θ, max_θ, θ_n)
    sliders, observables = getobservables(ranges...)
    all = lift(createobj, observables.morphz, photoreceptor_radius, absorption_coefficient)
    sc = Scene(scale_plot=false, limits=limits)
    vertices = lift(all) do all
        populateθ0!(all.layers)
        todo = [k for k in keys(all.layers) if k ∉ (:water, :mirror)]
        map(todo) do k
            getvertices(all.layers, k, nvertices)
        end
    end
    colors = lift(all) do all
        todo = [k for k in keys(all.layers) if k ∉ (:water, :mirror)]
        m, M = extrema(all.layers[k].tissue.ri for k in todo)
        Δ = M - m
        m -= 0.1Δ
        M += 0.1Δ
        if Δ == 0
            m = copy(M)
            M += 1
        end
        map(todo) do k
            layers2color(all.layers, k, m, M) 
        end
    end
    poly!(sc, vertices, color = colors)
    rays = lift(observables.distance, observables.aperture, all, observables.θ) do distance, aperture, all, θ
        l = Light(mm2μm(distance), aperture, all.rc, deg2rad(θ), all.nodalz)
        getrays(all.ous, l)
    end
    segments = lift(rays) do rays
        [ray[i].coord => ray[i+1].coord for ray in rays for i in 1:length(ray) - 1]
    end
    colors = lift(rays) do rays
        [mean((ray[i].color, ray[i+1].color)) for ray in rays for i in 1:length(ray) - 1]
    end
    linesegments!(sc, segments, color=colors)
    scatter!(sc, lift(a -> [Point2f0(0, a.nodalz)], all), markersize=10, color=:red)
    photoreceptors = lift(all) do all
        map([Scallops.RETINAS...]) do retina
            z1, z2 = Scallops.getzs(all.layers[retina], photoreceptor_radius)
            Rectangle{Float32}(-photoreceptor_radius, z1, 2photoreceptor_radius, z2 - z1)
        end
    end
    poly!(sc, photoreceptors, color=:green)
    axis = sc[Axis];
    axis[:names, :axisnames] = ("X (μm)", "Y (μm)");
    vbox(hbox(sliders...), sc)
end

