include(joinpath(dirname(pathof(Scallops)), "plotting.jl"))
using DataDeps, FileIO, Colors, CoordinateTransformations, ImageMorphology, Images, Interpolations
register(DataDep("scalloPupil", "a gif of scallop pupils", "https://vision-group-temporary.s3.amazonaws.com//scallop_pupil_Pm_green1_alt.gif", "4c8a194dd8a5df37833016d48a2f0f17cca321a8a6618984ddbc0a716b692063"))

folder = @datadep_str "scalloPupil"
file = first(readdir(folder))
gif = Gray.(load(joinpath(folder, file)))
bw = gif .< 0.1
bw = opening(bw)
lb = label_components(bw)
function fun(bb)
    XY = [zeros(2) for _ in bb]
    for (i,b) in enumerate(bb[2:2])
        xy, wh = b
        x,y = xy
        w,h = wh
        XY[i][1] = 2x + w
        XY[i][2] = 2y + h
    end
XY
end
wgif = [warp(gif[:,:,1], Translation(0,0))]
i = 1
bb = fun(component_boxes(lb[:,:,i]))
maxΔ = 0.0
for i in 2:size(lb, 3)
    global bb, maxΔ
    _bb = fun(component_boxes(lb[:,:,i]))
    Δ = -mean(bb .- _bb)
    maxΔ = max(maxΔ, abs.(Δ)...)
    tfm = Translation(Δ)
    push!(wgif, warp(gif[:,:,i], tfm))
end
w = cat(paddedviews(0, wgif...)..., dims = 3)
h = ceil(Int, maxΔ)
w = w[h+1:end-h, h+1:end-h, :]

bw = gif .> 0.1
bw = opening(bw)
lb = label_components(bw)

r = Vector{Float64}(undef, size(lb, 3))
for i in 1:size(lb,3)
    a = component_lengths(lb[:,:,i])[3]
    r[i] = sqrt.(a/π)
end
r .-= minimum(r)
r ./= maximum(r)
i = sortperm(r)
r = r[i]
sortedw = w[:,:,i]
itp = interpolate((1:size(w,1), 1:size(w,2), r), sortedw, Gridded(Linear()))
rl = range(0, 1, length=2size(w,3))
w2 = Array{Gray{Normed{UInt8,8}},3}(undef, size(w,1), size(w, 2), length(rl))
for i in 1:size(w,1), j in 1:size(w,2), (k,r) in enumerate(rl)
    w2[i,j,k] = itp(i,j,r)
end

min_aperture = 200
max_aperture = 400
aperture_n = size(w2,3)
min_morph = 0.8
max_morph = 1.48
morph_n = 100
min_distance = 0.1
max_distance = 1000
distance_n = 100
min_θ = -38
max_θ = 38
θ_n = max_θ - min_θ + 1
limits = FRect(-350,-210,700,420)
nvertices = 21
photoreceptor_radius = 5
absorption_coefficient = (dis_ret = 0.00667, pro_ret = 0.00667)
ranges = getranges(min_aperture, max_aperture, aperture_n, min_morph, max_morph, morph_n, min_distance, max_distance, distance_n, min_θ, max_θ, θ_n)
sliders, observables = getobservables(ranges...)
all2 = lift(createobj, observables.morphz, photoreceptor_radius, absorption_coefficient)
sc = Scene(scale_plot=false, limits=limits)
vertices2 = lift(all2) do all2
    populateθ0!(all2.layers)
    todo = [k for k in keys(all2.layers) if k ∉ (:water, :mirror)]
    map(todo) do k
        getvertices(all2.layers, k, nvertices)
    end
end
colors2 = lift(all2) do all2
    todo = [k for k in keys(all2.layers) if k ∉ (:water, :mirror)]
    m, M = extrema(all2.layers[k].tissue.ri for k in todo)
    Δ = M - m
    m -= 0.1Δ
    M += 0.1Δ
    if Δ == 0
        m = copy(M)
        M += 1
    end
    map(todo) do k
        layers2color(all2.layers, k, m, M) 
    end
end
poly!(sc, vertices2, color = colors2)
rays = lift(observables.distance, observables.aperture, all2, observables.θ) do distance, aperture, all2, θ
    l = Light(mm2μm(distance), aperture, all2.rc, deg2rad(θ), all2.nodalz)
    getrays(all2.ous, l)
end
segments = lift(rays) do rays
    [ray[i].coord => ray[i+1].coord for ray in rays for i in 1:length(ray) - 1]
end
colors2 = lift(rays) do rays
    [mean((ray[i].color, ray[i+1].color)) for ray in rays for i in 1:length(ray) - 1]
end
linesegments!(sc, segments, color=colors2)
axis = sc[Makie.Axis];
axis[:names, :axisnames] = ("X (μm)", "Y (μm)");
observables.morphz[] = 1.1
t = Node(1)
aspect_ratio = size(w2, 1)/size(w2, 2)
image!(sc, range(-400,400,length = size(w2,2)), range(250,250+800*aspect_ratio, length = size(w2,1)), lift(i -> w2[end:-1:1,:,i]', t), scale_plot=false, show_axis = false)
txtcolor = :black
rctcolor = :white
h = 15.0
y = 190.0
poly!(sc, FRect(-h,y-h,2h,2h), color = rctcolor)
text!(sc, "C", position = (0.0, y), align = (:center,  :center), textsize = 2h, color = txtcolor)
y = 60.0
poly!(sc, FRect(-h,y-h,2h,2h), color = rctcolor)
text!(sc, "L", position = (0.0, y), align = (:center,  :center), textsize = 2h, color = txtcolor)
y = -62.5
poly!(sc, FRect(-h,y-h,2h,2h), color = rctcolor)
text!(sc, "D", position = (0.0, y), align = (:center,  :center), textsize = 2h, color = txtcolor)
y = -150.0
poly!(sc, FRect(-h,y-h,2h,2h), color = rctcolor)
text!(sc, "P", position = (0.0, y), align = (:center,  :center), textsize = 2h, color = txtcolor)
y = -220.0
poly!(sc, FRect(-h,y-h,2h,2h), color = rctcolor)
text!(sc, "M", position = (0.0, y), align = (:center,  :center), textsize = 2h, color = txtcolor)
apertures = ranges.aperture#[end:-1:1]
record(sc, "pupil.mp4", 1:aperture_n) do i
    if i == 1
        observables.distance[] = maximum(ranges.distance)
    end
    observables.aperture[] = apertures[i]
    t[] = i
end

