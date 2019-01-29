using Scallops
using RayTraceEllipsoids
using DataFrames
using DataFramesMeta
using CSV

retina = :dis_ret
morphz = 1.1
photoreceptor_radius = 5
absorption_coefficient = (dis_ret = 0.00667, pro_ret = 0.00667)
distance = 1e5
layers, rc, ous, nodalz, focallengths = createobj(morphz, photoreceptor_radius, absorption_coefficient)
n = 5
apertures = range(200, 400, length=n)
df = DataFrame(method = String[], aperture = Float64[], Δθ = Float64[], fraction = Float64[], geometry = Float64[], sensitivity = Float64[])
for aperture in apertures
    Δθ = Scallops.getΔΘ(distance, aperture, rc, nodalz, ous, retina)
    cosa = cos(Δθ/2)
    fraction = Scallops.raytrace_fraction(distance, aperture, rc, ous, nodalz, cosa, retina)
    geometry = Scallops.unmodified_photons(distance, cosa, aperture, rc, nodalz)
    sensitivity = fraction*geometry
    push!(df, ("ray-traced", aperture, Δθ, fraction, geometry, sensitivity))
    Δθ = Scallops.getΔΘ(focallengths[retina], photoreceptor_radius)
    geometry = Scallops.modified_photons(aperture, Δθ)
    l = layers[retina].tissue.thick
    fraction = RayTraceEllipsoids.absorption(l, absorption_coefficient[retina])
    sensitivity = fraction*geometry
    push!(df, ("modified-Land", aperture, Δθ, fraction, geometry, sensitivity))
end

CSV.write("compare.csv", df)
for method in ("ray-traced", "modified-Land")
    d = @linq df |> where(:method .== method)
    CSV.write("$method.csv", d)
end

