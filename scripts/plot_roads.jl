using OpenCRG
using GLMakie
using Downloads: download

const REPO_ROOT = dirname(@__DIR__)
const REAL_ROAD_CACHE = joinpath(@__DIR__, "data", "real_roads")
const OUTDIR = joinpath(REPO_ROOT, "assets", "road_plots")
const REAL_ROAD_BASE_URL = "https://raw.githubusercontent.com/asam-ev/OpenCRG/master/crg-bin"

# Every plot uses true 1:1:1 real-world proportions -- nothing is scaled.
# Some roads' relief is genuinely tiny relative to their planar extent (a
# realistic road IS close to flat at true scale), which can read as sparse or
# make Z-axis ticks crowd together -- that's an honest property of the data,
# not something to hide by stretching the axes.
const ASPECT = :data

# (test fixture path relative to repo root, plot title)
const EXAMPLE_ROADS = [
    ("lib/OpenCRG/test/data/handmade_curved_minimalist.crg", "Handmade Curved Minimalist"),
    ("lib/OpenCRG/test/data/handmade_curved_banked_sloped.crg", "Handmade Curved, Banked & Sloped"),
    ("lib/OpenCRG/test/data/belgian_block.crg", "Belgian Block"),
    ("lib/OpenCRG/test/data/synthetic_hairpin.crg", "Synthetic Hairpin (180° turn)"),
    ("lib/OpenCRG/test/data/synthetic_end_anchored_2col.crg", "Synthetic End-Anchored"),
]

# (filename in asam-ev/OpenCRG's crg-bin/, plot title)
const REAL_ROADS = [
    ("country_road.crg", "Country Road (real-world survey, 3D Mapping Solutions GmbH)"),
    ("crg_refline_Hoki_HoeKi_Grafing.crg", "Hofolding–Höhenkirchen–Grafing Reference Line (real-world survey)"),
]

function ensure_real_road(filename)
    path = joinpath(REAL_ROAD_CACHE, filename)
    if !isfile(path)
        mkpath(REAL_ROAD_CACHE)
        url = "$REAL_ROAD_BASE_URL/$filename"
        println("Downloading $url ...")
        download(url, path)
    end
    return path
end

"""
Stride down a (nu, nv) grid so each dimension has at most `target` points.
Real-world surveyed roads can have hundreds of thousands of rows at cm
resolution -- far denser than a surface plot can usefully show or render
quickly, so we keep every `stride`-th sample rather than every sample.
"""
function downsample(X, Y, Z; target=400)
    nu, nv = size(Z)
    su = max(1, cld(nu, target))
    sv = max(1, cld(nv, target))
    return X[1:su:end, 1:sv:end], Y[1:su:end, 1:sv:end], Z[1:su:end, 1:sv:end]
end

function plot_road(path, title, outfile; target=400)
    println("== $title ==")
    println("  reading $path")
    data = read_crg(path)
    u, v, X, Y, Z = road_surface_grid(data)
    width = v[end] - v[1]
    println("  grid: $(size(Z)), width=$(width) m")
    fig = Figure(size=(1400, 1000))
    if width < 1.0
        # Sub-meter width means this is effectively a reference-line/path export
        # with no meaningful cross-section (e.g. a route-geometry-only survey) --
        # a "surface" over it would be a sub-pixel sliver, so plot the centerline
        # as a 3D curve instead, colored by distance traveled along the route.
        println("  width < 1m: not surface-plottable -- rendering as a 3D path instead")
        su = max(1, cld(length(u), target * target))
        idx = 1:su:length(u)
        centercol = clamp(cld(size(Z, 2), 2), 1, size(Z, 2))
        xs, ys, zs = X[idx, centercol], Y[idx, centercol], Z[idx, centercol]
        ax = Axis3(fig[1, 1]; xlabel="X (m)", ylabel="Y (m)", zlabel="Z (m)", aspect=ASPECT,
            title=title * "  [path — width = $(round(width, digits=2)) m, not a surface]")
        plt = GLMakie.lines!(ax, xs, ys, zs; color=u[idx], colormap=:viridis, linewidth=3)
        Colorbar(fig[1, 2], plt, label="Distance along route, u (m)")
    else
        println("  downsampling to <= $target per axis")
        Xd, Yd, Zd = downsample(X, Y, Z; target=target)
        println("  plotting $(size(Zd)) points")
        ax = Axis3(fig[1, 1]; xlabel="X (m)", ylabel="Y (m)", zlabel="Z (m)", title=title, aspect=ASPECT)
        plt = GLMakie.surface!(ax, Xd, Yd, Zd; colormap=:viridis)
        Colorbar(fig[1, 2], plt, label="Elevation (m)")
    end
    mkpath(dirname(outfile))
    save(outfile, fig)
    println("  saved $outfile")
end

for (relpath, title) in EXAMPLE_ROADS
    outfile = joinpath(OUTDIR, splitext(basename(relpath))[1] * ".png")
    plot_road(joinpath(REPO_ROOT, relpath), title, outfile)
end

for (filename, title) in REAL_ROADS
    path = ensure_real_road(filename)
    outfile = joinpath(OUTDIR, splitext(filename)[1] * ".png")
    plot_road(path, title, outfile; target=500)
end

println("Done. Plots written to $OUTDIR")
