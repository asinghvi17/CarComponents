# country_road_conceptcar_animation.jl — raytraced ControlledCountryRoadCarWithPowertrain video
# with the NVIDIA ConceptCar01 hero asset, rendered via OmniverseMakie (ovrtx RTX).
#
# Same model/solve as country_road_powertrain_animation.jl (whose helper functions this file
# reuses via include), but the visual REPLACES the MultibodyComponents primitive render with a
# `usdplot` of the ConceptCar USD driven by `bind_usd!`: per frame only 5 transforms are written
# (car root + 4 wheels: spin, front steer, suspension travel) instead of ~52 rebuilt meshes —
# ~18× faster than the MBC-render recording of the same solution (≈16 min for 50 s @ 24 fps).
# Scene extras: HDRI dome sky (visible background), procedurally textured asphalt + grass.
#
# Run (GPU box):
#   flock /tmp/omniversemakie-gpu.lock env DISPLAY=:0 JULIA_CUDA_USE_COMPAT=false \
#     OVRTX_LIBRARY_PATH=<...>/ovrtx/bin/libovrtx-dynamic.so \
#     julia +dyad-3.2.0-next.87 --project=. scripts/country_road_conceptcar_animation.jl
#
# Env knobs:
#   CONCEPTCAR_USD  path to ConceptCar01_Adjust.usd  (default: local digital-twins clone)
#   SKY_HDR         equirect HDR for the dome        (default: Kit packman table_mountain.hdr)
#   COUNTRY_RECORD  "0" to skip the video (stills + timing only; default records)
#   GRASS3D         "1" to add ~20k instanced 3D grass tufts (default off — texture-only ground)
#   SHOW_MBC        "1" to build the scene from MultibodyComponents.render so the multibody bits
#                   (suspension linkages, springs, axles, chassis rod) render inside/under the
#                   ConceptCar shell; the red MBC wheel cylinders are removed (the asset's wheels
#                   replace them at the same radius).  Slower: every kept MBC plot re-uploads its
#                   geometry each frame.  Writes to a separate _mbc video.
#   CAR_GLASS       "1" to swap the car's exterior shells (whole-car meshes exteriorDoorFlA +
#                   bonnet + carbon skirt + rubber seals — the asset has no per-panel meshes) to a
#                   ghost surface via a wrapper USDA that overrides their material:binding, so the
#                   suspension shows through the skin.  Pairs with SHOW_MBC=1.
#   GLASS_OPACITY   ghost-shell opacity (default 0.04; 0.0 = shells perfectly invisible)
#
# Assets/licensing: ConceptCar01 ships via Git LFS inside NVIDIA's digital-twins-for-fluid-
# simulation blueprint (internal dev/eval only — do not redistribute); the sky HDR ships in the
# Kit USD runtime resources (packman cache). Solve checkpoint + generated textures land in the
# gitignored scripts/data/.
#
# Load-bearing conventions (each cost an iteration to learn — see bug-reports + git history):
#   • The scene is Z-UP (OmniverseMakie stages author upAxis="Z", and the HDRI dome's poles
#     follow it): all sim data (Y-up) maps through p_scene = (x, −z, y).
#   • usdplot up = :z, NOT :y — the per-frame quaternion below is already the FULL asset→scene
#     rotation; up=:y would fold another +90°X in and bury the car nose-up in the road.
#   • Wheels spin at arc/r_RENDERED (the wheelbase-locked 0.677 scale makes the rendered wheel
#     r=0.26 m vs the sim's 0.20 m; sim wheel speed would look like slip).
#   • Road-ribbon UVs use the center-line ARC LENGTH (x would stretch the texture on curves).
#   • The CRG lane is ASYMMETRIC: v ∈ [−1.57, +1.07] m — grass edges are per-side.

using Serialization
using LinearAlgebra
import OmniverseMakie as OM

include(joinpath(@__DIR__, "country_road_powertrain_animation.jl"))   # model + road + helpers
const MK = OM.Makie

const ROOT = dirname(@__DIR__)
const DATA = joinpath(@__DIR__, "data")                                # gitignored
mkpath(DATA)
const CAR = get(ENV, "CONCEPTCAR_USD",
    "/home/juliahub/temp/digital-twins-for-fluid-simulation/stages/layers/CarHero/ConceptCar/ConceptCar01_Adjust.usd")
const SKY = get(ENV, "SKY_HDR",
    "/home/juliahub/.cache/packman/chk/usd.py312.manylinux_2_35_x86_64.stock.release/0.25.02.kit.8-gl.16788+c1c423f2/resources/Lights/table_mountain.hdr")
const SOL_CACHE = joinpath(DATA, "country_conceptcar_sol.jls")
const SHOW_MBC  = get(ENV, "SHOW_MBC", "0") == "1"
const CAR_GLASS = get(ENV, "CAR_GLASS", "0") == "1"
# Ghost-shell opacity for CAR_GLASS: 0.0 = the swapped shells are perfectly invisible; ~0.03–0.05
# leaves a faint hint of the body shape.  (Fractional OmniPBR opacity, not refractive glass —
# refraction/Fresnel read as "smoked glass" and hide the suspension.)
const GLASS_OPACITY = parse(Float64, get(ENV, "GLASS_OPACITY", "0.04"))
const VIDEO_OUT = joinpath(ROOT, "assets",
    "country_road_conceptcar" * (SHOW_MBC ? "_mbc" : "") * (CAR_GLASS ? "_glass" : "") * "_rt.mp4")
logmsg(m) = (println("[", round(time(); digits = 1), "] ", m); flush(stdout))

isfile(CAR) || error("ConceptCar asset missing: $CAR (set CONCEPTCAR_USD)")
isfile(SKY) || error("sky HDR missing: $SKY (set SKY_HDR)")

# ================================= solve (or load cache) =========================================
# The symbolic initialization-problem construction hangs on some machines (SymbolicUtils
# hashconsing blowup on controller-feedback models — bug-reports/init-hashconsing-hang/), so the
# problem is built with build_initializeprob=false, the algebraic unknowns are seeded, and
# BrownFullBasicInit() initializes numerically at solve time.
function algebraic_seeds(sys, start_z)
    corners = [sys.excited_suspension_fr, sys.excited_suspension_fl,
               sys.excited_suspension_br, sys.excited_suspension_bl]
    seeds = Pair[]
    for c in corners
        for (path, val) in [((:wheel, :wheeljoint, :delta_0), nothing),   # vector: [1],[2],[3] => 0
                            ((:rotational_losses, :w_rel), 0.0),
                            ((:rotational_losses, :d), 0.01),
                            ((:suspension, :r2, :w), 0.0),
                            ((:wheel_lateral_position,), start_z)]
            obj = c
            ok = true
            for p in path
                try
                    obj = getproperty(obj, p)
                catch
                    ok = false
                    break
                end
            end
            ok || continue
            if val === nothing
                push!(seeds, obj[1] => 0.0, obj[2] => 0.0, obj[3] => 0.0)
            else
                push!(seeds, obj => val)
            end
        end
    end
    return seeds
end

function simulate_noinit(sys, road_surface, center_z_profile, center_heading_profile;
        speed, target_speed, tstop, save_dt, baseline = nothing, start_x = 2.0)
    defs = country_initial_conditions(sys, road_surface, center_z_profile, center_heading_profile;
        speed = speed, baseline = baseline, start_x = start_x)
    append!(defs, [sys.target_speed => target_speed, sys.controller.target_speed => target_speed])
    append!(defs, algebraic_seeds(sys, center_z_profile(start_x)))
    prob = ODEProblem(sys, defs, (0.0, tstop);
        optimize = :basic, saveat = save_dt,
        warn_initialize_determined = false, build_initializeprob = false)
    sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff());
        initializealg = BrownFullBasicInit(), abstol = 1e-6, reltol = 1e-6, dtmax = 0.02)
    SciMLBase.successful_retcode(sol) || error("solve failed: $(sol.retcode), t_end=$(sol.t[end])")
    return sol
end

function solve_or_load()
    if isfile(SOL_CACHE)
        logmsg("loading solve cache $SOL_CACHE ...")
        ck = deserialize(SOL_CACHE)
        return ck.sol, ck.road_data
    end
    logmsg("preparing country-road data...")
    road_data = CarComponents.country_road_curved_data(longitudinal_spacing = 0.5, lateral_stride = 8)
    road_surface = CarComponents.opencrg_curved_surface_interpolator(road_data)
    center_z_profile = CarComponents.opencrg_curved_center_z_interpolator(road_data)
    center_heading_profile = CarComponents.opencrg_curved_heading_interpolator(road_data)
    model, sys = compile_controlled_country_model(road_surface, center_z_profile, center_heading_profile;
        target_speed = 6.0)
    logmsg("compiled: $(length(unknowns(sys))) unknowns; settling static baseline...")
    bsol = simulate_noinit(sys, road_surface, center_z_profile, center_heading_profile;
        speed = 0.0, target_speed = 0.0, tstop = 4.0, save_dt = 1 / 60)
    tb = bsol.t[end]
    body_x = bsol(tb; idxs = sys.back_front.body.r_0[1])
    body_z = bsol(tb; idxs = sys.back_front.body.r_0[3])
    baseline = (
        body_y = bsol(tb; idxs = sys.back_front.body.r_0[2]),
        body_ay = bsol(tb; idxs = sys.back_front.body.a_0[2]),
        body_x = body_x, body_z = body_z,
        road_body_center = road_surface(body_x, body_z),
        fr_s = bsol(tb; idxs = sys.excited_suspension_fr.suspension.springdamper.s),
        fl_s = bsol(tb; idxs = sys.excited_suspension_fl.suspension.springdamper.s),
        br_s = bsol(tb; idxs = sys.excited_suspension_br.suspension.springdamper.s),
        bl_s = bsol(tb; idxs = sys.excited_suspension_bl.suspension.springdamper.s),
    )
    logmsg("moving run (target 6 m/s, 50 s)...")
    sol = simulate_noinit(sys, road_surface, center_z_profile, center_heading_profile;
        speed = 0.0, target_speed = 6.0, tstop = 50.0, save_dt = 1 / 24, baseline = baseline)
    serialize(SOL_CACHE, (; sol, baseline, road_data))
    logmsg("solve cached: $SOL_CACHE ($(round(filesize(SOL_CACHE) / 1e6; digits = 1)) MB)")
    return sol, road_data
end

sol, road_data = solve_or_load()
road_surface = CarComponents.opencrg_curved_surface_interpolator(road_data)
center_z_profile = CarComponents.opencrg_curved_center_z_interpolator(road_data)
center_heading_profile = CarComponents.opencrg_curved_heading_interpolator(road_data)
model, sys = compile_controlled_country_model(road_surface, center_z_profile, center_heading_profile;
    target_speed = 6.0)
logmsg("model ready: $(length(sol.t)) samples, t=[$(sol.t[1]), $(sol.t[end])]")

# ============================ per-frame signals (extracted once) =================================
series(var) = collect(sol[var])
fa = (x = series(sys.back_front.frame_a.r_0[1]),          # front-axle center (anchor)
      y = series(sys.back_front.frame_a.r_0[2]),
      z = series(sys.back_front.frame_a.r_0[3]))
fb = (x = series(sys.back_front.frame_b.r_0[1]),          # rear-axle center
      y = series(sys.back_front.frame_b.r_0[2]),
      z = series(sys.back_front.frame_b.r_0[3]))
axR = (x = series(sys.front_axle.frame_a.r_0[1]),         # front axle RIGHT end (fr)
       y = series(sys.front_axle.frame_a.r_0[2]),
       z = series(sys.front_axle.frame_a.r_0[3]))
axL = (x = series(sys.front_axle.frame_b.r_0[1]),         # front axle LEFT end (fl)
       y = series(sys.front_axle.frame_b.r_0[2]),
       z = series(sys.front_axle.frame_b.r_0[3]))
steer = series(sys.controller.steer_angle)
s_fr = series(sys.excited_suspension_fr.suspension.springdamper.s)
s_fl = series(sys.excited_suspension_fl.suspension.springdamper.s)
s_br = series(sys.excited_suspension_br.suspension.springdamper.s)
s_bl = series(sys.excited_suspension_bl.suspension.springdamper.s)
N = length(sol.t)

arc = zeros(N)                                            # distance rolled, for wheel spin
for i in 2:N
    arc[i] = arc[i-1] + hypot(fa.x[i] - fa.x[i-1], fa.z[i] - fa.z[i-1])
end

# steer-sign calibration: +R_y in the wheel-parent frame steers LEFT; left turn about sim +y has
# (fwd_i × fwd_{i+1})·ŷ > 0 — calibrate rather than hardcode the sim's sign convention.
fwdx = fa.x .- fb.x; fwdz = fa.z .- fb.z
lrate = zeros(N)
for i in 1:N-1
    lrate[i] = fwdz[i] * fwdx[i+1] - fwdx[i] * fwdz[i+1]
end
const K_STEER = sum(steer .* lrate) >= 0 ? 1.0f0 : -1.0f0
logmsg("steer sign: K_STEER=$(K_STEER)")

# ====================== ConceptCar01 geometry constants (all measured) ===========================
# Source: bug-reports/car-mesh-sourcing/vehicle_dimensions.json (pxr measurement of the asset).
# Asset: Y-up, authored cm, nose −z, origin ~front axle; wheelbase 2.9538 m, mean track 1.862 m.
# The model's COUNTRY_ROAD_TARGET_* constants are this geometry scaled to the sim's 1.48 m
# contact track, so deriving the render scale from axle_spacing makes wheelbase AND track exact
# (rendered track = 1.862 × 0.7949 = 1.480 m); wheel radii keep the asset's front/rear stagger
# (0.304/0.310 m vs the sim's single 0.3076 — sub-1% and invisible).
const S_SCALE = COUNTRY_ROAD_TARGET_AXLE_SPACING / 2.9538   # ≈ 0.79494
const S_TOT   = S_SCALE * 0.01         # asset cm → world m
const R_VIS_F = 0.3827 * S_SCALE       # rendered wheel radii [m]
const R_VIS_B = 0.3900 * S_SCALE
const C_FWD   = 0.0936 * S_SCALE       # asset origin sits this far AHEAD of the front axle
δ = fa.y[1] - road_surface(fa.x[1], fa.z[1])   # settled front-axle height above road
logmsg("ride height: δ=$(round(δ; digits = 4)) m")

function body_pose(i)                                     # sim-frame R (columns fwd/up/right), O
    X = normalize([fwdx[i], fa.y[i] - fb.y[i], fwdz[i]])
    lat = [axR.x[i] - axL.x[i], axR.y[i] - axL.y[i], axR.z[i] - axL.z[i]]
    Z = normalize(lat .- (dot(lat, X) .* X))
    Y = cross(Z, X)
    R = [X Y Z]
    O = [fa.x[i], fa.y[i], fa.z[i]] .+ R * [C_FWD, -δ, 0.0]
    return R, O
end

function rot_to_quat(R)                                   # self-checked vs Makie's convention
    t = R[1,1] + R[2,2] + R[3,3]
    if t > 0
        s = sqrt(t + 1.0) * 2
        w = 0.25s; x = (R[3,2]-R[2,3])/s; y = (R[1,3]-R[3,1])/s; z = (R[2,1]-R[1,2])/s
    elseif R[1,1] > R[2,2] && R[1,1] > R[3,3]
        s = sqrt(1.0 + R[1,1] - R[2,2] - R[3,3]) * 2
        w = (R[3,2]-R[2,3])/s; x = 0.25s; y = (R[1,2]+R[2,1])/s; z = (R[1,3]+R[3,1])/s
    elseif R[2,2] > R[3,3]
        s = sqrt(1.0 + R[2,2] - R[1,1] - R[3,3]) * 2
        w = (R[1,3]-R[3,1])/s; x = (R[1,2]+R[2,1])/s; y = 0.25s; z = (R[2,3]+R[3,2])/s
    else
        s = sqrt(1.0 + R[3,3] - R[1,1] - R[2,2]) * 2
        w = (R[2,1]-R[1,2])/s; x = (R[1,3]+R[3,1])/s; y = (R[2,3]+R[3,2])/s; z = 0.25s
    end
    q = MK.Quaternionf(x, y, z, w)
    err = maximum(abs.(Matrix(MK.rotationmatrix4(q))[1:3, 1:3] .- R))
    err < 1e-3 || error("quaternion convention mismatch (err=$err)")
    return q
end

# sim (Y-up) → scene (Z-up): p_scene = (x, −z, y) = R_x(+90°)·p_sim (det +1, no mirror)
s3(x, y, z) = MK.Vec3f(x, -z, y)
s3(v) = s3(v[1], v[2], v[3])
const Q_X90  = MK.qrotation(MK.Vec3f(1, 0, 0), Float32(π / 2))
const Q_BASE = MK.qrotation(MK.Vec3f(0, 1, 0), Float32(-π / 2))   # asset nose(−z) → sim +x

# ====================== procedural textures (deterministic, cached in DATA) ======================
function _noise_octave(nextf, px, cells)                  # periodic bilinear value noise
    g = [nextf() for i in 1:cells, j in 1:cells]
    img = Matrix{Float32}(undef, px, px)
    cs = px / cells
    for j in 1:px, i in 1:px
        u = (i - 1) / cs; v = (j - 1) / cs
        i0 = floor(Int, u); j0 = floor(Int, v)
        fu = Float32(u - i0); fv = Float32(v - j0)
        su = fu * fu * (3 - 2fu); sv = fv * fv * (3 - 2fv)
        a = g[i0 % cells + 1, j0 % cells + 1]
        b = g[(i0 + 1) % cells + 1, j0 % cells + 1]
        c = g[i0 % cells + 1, (j0 + 1) % cells + 1]
        d = g[(i0 + 1) % cells + 1, (j0 + 1) % cells + 1]
        img[i, j] = (a * (1 - su) + b * su) * (1 - sv) + (c * (1 - su) + d * su) * sv
    end
    return img
end
_lcg(seed) = let state = Ref(seed)
    () -> (state[] = state[] * 0x5851f42d4c957f2d + 0x14057b7ef767814f;
           Float32((state[] >> 40) % UInt64(2)^24) / 2.0f0^24)
end

const GRASS_PNG = joinpath(DATA, "grass_albedo.png")
function make_grass!(path; px = 1024)
    nextf = _lcg(UInt64(0x9e3779b97f4a7c15))
    n  = 0.50f0 .* _noise_octave(nextf, px, 128) .+ 0.30f0 .* _noise_octave(nextf, px, 32) .+
         0.20f0 .* _noise_octave(nextf, px, 8)
    lo = _noise_octave(nextf, px, 4)                      # broad dry patches
    img = map(n, lo) do t, d
        g1 = (0.16f0, 0.25f0, 0.09f0); g2 = (0.34f0, 0.44f0, 0.16f0); dry = (0.38f0, 0.34f0, 0.18f0)
        r = g1[1] + (g2[1] - g1[1]) * t; g = g1[2] + (g2[2] - g1[2]) * t; b = g1[3] + (g2[3] - g1[3]) * t
        w = clamp((d - 0.62f0) * 2.2f0, 0f0, 0.28f0)
        MK.RGBf(r + (dry[1] - r) * w, g + (dry[2] - g) * w, b + (dry[3] - b) * w)
    end
    OM.PNGFiles.save(path, img)
end

const ASPHALT_PNG = joinpath(DATA, "asphalt_albedo.png")
function make_asphalt!(path; px = 1024)
    nextf = _lcg(UInt64(0x51e2d24c3a7f9b11))
    patch = 0.6f0 .* _noise_octave(nextf, px, 64) .+ 0.4f0 .* _noise_octave(nextf, px, 16)
    img = Matrix{MK.RGBf}(undef, px, px)
    for j in 1:px, i in 1:px
        g = 0.135f0 + 0.05f0 * (patch[i, j] - 0.5f0)      # binder tone
        g *= 0.75f0 + 0.5f0 * nextf()                     # per-pixel aggregate speckle
        nextf() < 0.004f0 && (g += 0.20f0 + 0.15f0 * nextf())   # sparse quartz glints
        g = clamp(g, 0.0f0, 1.0f0)
        img[i, j] = MK.RGBf(min(g * 1.04f0, 1), g, g * 0.96f0)
    end
    OM.PNGFiles.save(path, img)
end

isfile(GRASS_PNG) || make_grass!(GRASS_PNG)
isfile(ASPHALT_PNG) || make_asphalt!(ASPHALT_PNG)

# ============================== scene (Z-up) =====================================================
OM.activate!(warmup = 24, accumulate_across_frames = true, background = :domelight)
lights = MK.AbstractLight[
    MK.DirectionalLight(MK.RGBf(1.8, 1.75, 1.6), s3(-0.45, -1.0, -0.3)),   # sun (dome adds ambient)
]
# W3: SIM-frame point → authored point.  The default scene authors the Z-up world directly (s3).
# The SHOW_MBC scene keeps content in SIM (Y-up) coordinates — that is what the MBC plots are
# authored in — and carries a single Q_X90 root rotation on the scene instead (plots inherit the
# scene transformation), so the dome's Z-up poles still line up.  The camera lives in world
# space either way, so its s3 mapping below is branch-independent.
if SHOW_MBC
    fig_mbc, time_obs, ls = MultibodyComponents.render(model, sol, sol.t[1];
        slider = false, show_axis = false, size = (1280, 720), lights = lights)
    scene = ls.scene
    MK.rotate!(scene, Q_X90)
    axplots = filter(p -> p isa MK.Axis3D, copy(scene.plots))
    foreach(p -> delete!(scene, p), axplots)
    # drop the red plots: the SlippingWheel cylinders (color default [1,0,0,1] — the ConceptCar
    # wheels replace them at the same radius) and the world x-axis arrow
    redplots = filter(copy(scene.plots)) do p
        c = try p.color[] catch; nothing end
        c isa MK.Colors.Colorant && Float32(MK.Colors.red(c)) > 0.9f0 &&
            Float32(MK.Colors.green(c)) < 0.15f0 && Float32(MK.Colors.blue(c)) < 0.15f0
    end
    foreach(p -> delete!(scene, p), redplots)
    # The dyad model colors the chassis/axle rods 30%-alpha gray — invisible through the glass
    # body.  Force every kept MBC plot bright and OPAQUE so the hardware reads as a skeleton.
    recolored = Ref(0)                                    # Ref: top-level for loops are soft scope
    for p in copy(scene.plots)
        c = try p.color[] catch; nothing end
        c isa MK.Colors.Colorant || continue
        r = Float32(MK.Colors.red(c)); g = Float32(MK.Colors.green(c)); b = Float32(MK.Colors.blue(c))
        p.color[] = MK.RGBAf(min(1, 0.3f0 + r), min(1, 0.3f0 + g), min(1, 0.3f0 + b), 1.0f0)
        recolored[] += 1
    end
    logmsg("MBC scene: removed $(length(axplots)) axis + $(length(redplots)) red plots; " *
           "$(length(scene.plots)) kept ($(recolored[]) recolored opaque)")
    W3(x, y, z) = MK.Point3f(x, y, z)
else
    scene = MK.Scene(size = (1280, 720); lights = lights)
    MK.cam3d!(scene; center = false)
    W3(x, y, z) = MK.Point3f(x, -z, y)
end

# Road ribbon from the SIM'S OWN interpolator functions: the 2D z(x) fits diverge from the CRG
# grid's v=0 column by meters on curves, and the car drives the FUNCTIONS' road.
vlo, vhi = extrema(road_data.v_axis)
xs = collect(range(0.0, 285.0; step = 0.5))
vs = collect(range(vlo, vhi; length = 9))
Xm = [xs[i] - vs[j] * sin(center_heading_profile(xs[i])) for i in eachindex(xs), j in eachindex(vs)]
Zm = [center_z_profile(xs[i]) + vs[j] * cos(center_heading_profile(xs[i])) for i in eachindex(xs), j in eachindex(vs)]
Ym = [road_surface(Xm[i, j], Zm[i, j]) for i in eachindex(xs), j in eachindex(vs)]

# per-vertex UVs (same vertex order as road_mesh: j outer, i inner); u = arc length / tile.
# Arrays are SIM-frame; W3 maps each vertex into the authored scene.
function road_mesh_uv(Xm, Ym, Zm, U, V)
    points = MK.Point3f[]
    uv = MK.Vec2f[]
    for j in axes(Xm, 2), i in axes(Xm, 1)
        push!(points, W3(Xm[i, j], Ym[i, j], Zm[i, j]))
        push!(uv, MK.Vec2f(U[i], V[j]))
    end
    nx = size(Xm, 1)
    faces = GLMakie.GeometryBasics.GLTriangleFace[]
    for j in 1:(size(Xm, 2) - 1), i in 1:(nx - 1)
        a = (j - 1) * nx + i; b = a + 1; c = j * nx + i; d = c + 1
        push!(faces, GLMakie.GeometryBasics.GLTriangleFace(a, c, b))
        push!(faces, GLMakie.GeometryBasics.GLTriangleFace(b, c, d))
    end
    return GLMakie.GeometryBasics.Mesh(points, faces; uv = uv)
end

const ASPHALT_TILE = 2.0
czl = [center_z_profile(x) for x in xs]
s_arc = zeros(length(xs))
for i in 2:length(xs)
    s_arc[i] = s_arc[i-1] + hypot(xs[i] - xs[i-1], czl[i] - czl[i-1])
end
U = Float32.(s_arc ./ ASPHALT_TILE)
V = Float32.((vs .- vlo) ./ ASPHALT_TILE)
MK.mesh!(scene, road_mesh_uv(Xm, Ym, Zm, U, V);
         color = :white, material = (; base_color_texture = ASPHALT_PNG,
                                       roughness = 0.9f0, metallic = 0.0f0))

# grass ground: one textured quad just under the road, UVs tiled every ~5 m
hg = minimum(Ym) - 0.04                                   # plane height, sim frame
xlo, xhi = extrema(Xm); zlo, zhi = extrema(Zm)
gx0, gx1 = xlo - 120, xhi + 120
gz0, gz1 = zlo - 120, zhi + 120
ntx = round((gx1 - gx0) / 5); nty = round((gz1 - gz0) / 5)
gpts = [W3(gx0, hg, gz0), W3(gx1, hg, gz0), W3(gx1, hg, gz1), W3(gx0, hg, gz1)]
guv  = MK.Vec2f[(0, 0), (ntx, 0), (ntx, nty), (0, nty)]
gnrm = [SHOW_MBC ? MK.Vec3f(0, 1, 0) : MK.Vec3f(0, 0, 1) for _ in 1:4]
gfcs = GLMakie.GeometryBasics.GLTriangleFace[(1, 2, 3), (1, 3, 4)]
MK.mesh!(scene, GLMakie.GeometryBasics.Mesh(gpts, gfcs; uv = guv, normal = gnrm); color = :white,
         material = (; base_color_texture = GRASS_PNG, roughness = 0.95f0, metallic = 0.0f0))

# optional 3D grass verge (~20k instanced tufts, one PointInstancer, zero per-frame cost but
# ~35% slower RTX steps).  Off by default — the texture-only ground was the preferred look.
# (Built in Z-up scene coordinates, so not available in the SHOW_MBC sim-frame scene.)
if get(ENV, "GRASS3D", "0") == "1" && SHOW_MBC
    @warn "GRASS3D=1 is not supported with SHOW_MBC=1 — skipping the 3D grass."
elseif get(ENV, "GRASS3D", "0") == "1"
    grand = _lcg(UInt64(0x243f6a8885a308d3))
    function tuft_marker(; nblades = 6)
        pts = MK.Point3f[]; nrm = MK.Vec3f[]; fcs = GLMakie.GeometryBasics.GLTriangleFace[]
        for b in 1:nblades
            θ = 2π * (b - 1) / nblades + 0.9f0 * (grand() - 0.5f0)
            bend = 0.15f0 + 0.35f0 * grand()
            dir = MK.Vec3f(cos(θ), sin(θ), 0); side = MK.Vec3f(-sin(θ), cos(θ), 0)
            ts = (0.0f0, 0.55f0, 1.0f0); ws = (0.05f0, 0.03f0, 0.004f0)
            st = [MK.Point3f((bend * t^2) * dir + MK.Vec3f(0, 0, t * (1 - 0.3f0 * bend * t^2))) for t in ts]
            tang = [normalize(st[min(k + 1, 3)] - st[max(k - 1, 1)]) for k in 1:3]
            n = [normalize(cross(side, tang[k])) for k in 1:3]
            for flip in (1.0f0, -1.0f0)                   # both windings, SEPARATED ±1.5 mm:
                i0 = length(pts)                          # coplanar twins z-fight → black blades
                for k in 1:3
                    off = 0.0015f0 * flip * n[k]
                    push!(pts, st[k] - ws[k] * side + off); push!(nrm, flip * n[k])
                    push!(pts, st[k] + ws[k] * side + off); push!(nrm, flip * n[k])
                end
                for k in (0, 2)
                    a, b2, c, d = i0 + k + 1, i0 + k + 2, i0 + k + 3, i0 + k + 4
                    if flip > 0
                        push!(fcs, GLMakie.GeometryBasics.GLTriangleFace(a, b2, c))
                        push!(fcs, GLMakie.GeometryBasics.GLTriangleFace(b2, d, c))
                    else
                        push!(fcs, GLMakie.GeometryBasics.GLTriangleFace(a, c, b2))
                        push!(fcs, GLMakie.GeometryBasics.GLTriangleFace(b2, c, d))
                    end
                end
            end
        end
        return GLMakie.GeometryBasics.Mesh(pts, fcs; normal = nrm)
    end
    # density bands (band start, band end, tufts/m²) outward from each road edge
    bands = ((0.05, 3.0, 5.5), (3.0, 7.0, 2.2), (7.0, 14.0, 0.8), (14.0, 25.0, 0.18))
    gpos = MK.Point3f[]; gsiz = MK.Vec3f[]; grot = MK.Quaternionf[]; gcol = MK.RGBf[]
    for i in 1:(length(xs) - 1)
        seg = hypot(xs[i+1] - xs[i], czl[i+1] - czl[i])
        hd = center_heading_profile(xs[i])
        edge_hL = road_surface(xs[i] - (vhi - 0.05) * sin(hd), czl[i] + (vhi - 0.05) * cos(hd))
        edge_hR = road_surface(xs[i] - (vlo + 0.05) * sin(hd), czl[i] + (vlo + 0.05) * cos(hd))
        for (d0, d1, dens) in bands, sgn in (1.0, -1.0)
            μ = dens * (d1 - d0) * seg
            cnt = floor(Int, μ) + (grand() < μ - floor(μ) ? 1 : 0)
            for _ in 1:cnt
                d = d0 + (d1 - d0) * grand()
                v = sgn > 0 ? vhi + d : vlo - d           # v_axis is ASYMMETRIC (−1.57…+1.07):
                xg = xs[i] + seg * grand() - v * sin(hd)  # per-side edges keep grass off the lane
                zg2 = czl[i] + v * cos(hd)
                eh = sgn > 0 ? edge_hL : edge_hR
                h = eh + (hg - eh) * clamp(d / 6, 0, 1)   # blend edge height → plane height
                push!(gpos, MK.Point3f(xg, -zg2, h - 0.03))
                ht = 0.16f0 + 0.18f0 * grand()
                push!(gsiz, MK.Vec3f(ht * (0.8f0 + 0.35f0 * grand()), ht * (0.8f0 + 0.35f0 * grand()), ht))
                push!(grot, MK.qrotation(MK.Vec3f(0, 0, 1), 2f0π * grand()))
                t = grand()
                c = MK.RGBf(0.16f0 + 0.20f0 * t, 0.26f0 + 0.21f0 * t, 0.08f0 + 0.08f0 * t)
                if grand() < 0.15f0
                    w = 0.5f0 * grand()
                    c = MK.RGBf(c.r + (0.42f0 - c.r) * w, c.g + (0.38f0 - c.g) * w, c.b + (0.20f0 - c.b) * w)
                end
                push!(gcol, c)
            end
        end
    end
    MK.meshscatter!(scene, gpos; marker = tuft_marker(), markersize = gsiz, rotation = grot,
                    color = gcol)
    logmsg("3D grass: $(length(gpos)) tufts")
end

# center reference line + driven path, thin WORLD-SPACE lines (OmniverseMakie linewidth = meters)
MK.lines!(scene, [W3(xs[i], road_surface(xs[i], czl[i]) + 0.03, czl[i]) for i in eachindex(xs)];
          color = :deepskyblue, linewidth = 0.06)
MK.lines!(scene, [W3(fa.x[i], road_surface(fa.x[i], fa.z[i]) + 0.05, fa.z[i]) for i in 1:N];
          color = :orange, linewidth = 0.05)

# ================================== the car ======================================================
# CAR_GLASS: a generated wrapper layer references the car and rebinds the two exterior paint
# shells to a self-contained thin-walled OmniGlass (resolved from ovrtx's core MDL library).
# The wrapper's defaultPrim reproduces the car's /World, so bind targets are unchanged.
CAR_EFFECTIVE = CAR
if CAR_GLASS
    wrap = joinpath(DATA, "conceptcar_glass_wrapper.usda")
    write(wrap, """
#usda 1.0
(
    defaultPrim = "World"
)

def Xform "World" (
    prepend references = @$(CAR)@
)
{
    def Scope "GlassSwap"
    {
        def Material "BodyGlass"
        {
            token outputs:mdl:surface.connect = </World/GlassSwap/BodyGlass/Shader.outputs:out>

            def Shader "Shader"
            {
                uniform token info:implementationSource = "sourceAsset"
                uniform asset info:mdl:sourceAsset = @OmniPBR.mdl@
                uniform token info:mdl:sourceAsset:subIdentifier = "OmniPBR"
                color3f inputs:diffuse_color_constant = (0.45, 0.75, 0.55)
                float inputs:reflection_roughness_constant = 1.0
                float inputs:metallic_constant = 0.0
                bool inputs:enable_opacity = 1
                float inputs:opacity_constant = $(GLASS_OPACITY)
                token outputs:out
            }
        }
    }

    over "root"
    {
        over "Body"
        {
            rel material:binding = </World/GlassSwap/BodyGlass> (
                bindMaterialAs = "strongerThanDescendants"
            )
        }
    }
}
""")
    global CAR_EFFECTIVE = wrap
    logmsg("glass wrapper written: $wrap")
end
car = OM.usdplot!(scene, CAR_EFFECTIVE; up = :z,
    bbox = MK.Rect3f(MK.Point3f(-111.4, -0.1, -100.7), MK.Vec3f(222.7, 137.4, 505.0)))
MK.scale!(car, S_TOT, S_TOT, S_TOT)

# Wheel groups /root/Wheels/wheel_{DF,PF,DB,PB} (D=driver=asset−x=fl/bl, P=passenger=fr/br).
# Their parent frame is cm-scale, x=LEFT, y=up, z=nose (right-handed, verified from the composed
# root·Wheels matrices): +R_x = forward roll, +R_y = steer left.  omni:xform REPLACES the local
# matrix, so the authored local L (P-side wheels are mirrored, det −1) is re-applied inside, and
# spin/steer compose about the measured hub (prim origins are NOT the hubs).
L_D = Float32[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
L_P = Float32[-0.999995 0 0.003051 0; 0 1 0 0; 0.003051 0 0.999995 0; 0 0 0 1]
wheel_mat(hub, dy, st, sp, L) =
    MK.translationmatrix(hub + MK.Vec3f(0, dy, 0)) * MK.rotationmatrix_y(st) *
    MK.rotationmatrix_x(sp) * MK.translationmatrix(-hub) * L
WHEELS = [   # (bind target, hub in parent cm, steerable, L, rendered radius m, suspension series)
    ("/root/Wheels/wheel_DF", MK.Vec3f(91.68, 38.37, 141.04), true,  L_D, R_VIS_F, s_fl),
    ("/root/Wheels/wheel_PF", MK.Vec3f(-91.7, 38.4, 141.2),   true,  L_P, R_VIS_F, s_fr),
    ("/root/Wheels/wheel_DB", MK.Vec3f(94.5, 38.9, -154.2),   false, L_D, R_VIS_B, s_bl),
    ("/root/Wheels/wheel_PB", MK.Vec3f(-94.4, 38.9, -154.2),  false, L_P, R_VIS_B, s_br),
]
wobs = Dict(tgt => OM.bind_usd!(car, tgt, wheel_mat(hub, 0.0f0, 0.0f0, 0.0f0, L))
            for (tgt, hub, _, L, _, _) in WHEELS)

function set_frame!(i)
    R, O = body_pose(i)
    # authored car rotation: the MBC scene already carries Q_X90 at the root, the default scene
    # needs it composed in here — the WORLD-space result is identical either way (the assertion
    # below checks that world result)
    MK.rotate!(car, SHOW_MBC ? rot_to_quat(R) * Q_BASE : Q_X90 * rot_to_quat(R) * Q_BASE)
    MK.translate!(car, MK.Vec3f(W3(O...)))
    SHOW_MBC && (time_obs[] = sol.t[i])                   # drive the MBC plots
    for (tgt, hub, steerable, L, rvis, sser) in WHEELS
        dy = Float32(-(sser[i] - sser[1]) / S_TOT)        # spring shortens → wheel UP (parent cm)
        st = steerable ? K_STEER * Float32(steer[i]) : 0.0f0
        sp = Float32(mod2pi(arc[i] / rvis))               # +R_x = forward roll
        wobs[tgt][] = wheel_mat(hub, dy, st, sp, L)
    end
    h = atan(fwdz[i], fwdx[i])
    mx, mz = (fa.x[i] + fb.x[i]) / 2, (fa.z[i] + fb.z[i]) / 2
    rh = road_surface(mx, mz)
    eye    = (mx - 6.2cos(h) - 2.3sin(h), rh + 2.5, mz - 6.2sin(h) + 2.3cos(h))
    lookat = (mx + 3.0cos(h), rh + 0.55, mz + 3.0sin(h))
    MK.update_cam!(scene, s3(eye...), s3(lookat...), MK.Vec3f(0, 0, 1))
    return nothing
end
set_frame!(1)

# fail-fast: the composed rotation must map asset nose(−z)→scene fwd and asset up(y)→scene +z
let (R1, _) = body_pose(1)
    M = Matrix(MK.rotationmatrix4(Q_X90 * rot_to_quat(R1) * Q_BASE))[1:3, 1:3]
    h1 = atan(fwdz[1], fwdx[1])
    fwd_scene = [cos(h1), -sin(h1), 0.0]
    isapprox(M * [0, 0, -1], fwd_scene; atol = 0.05) ||
        error("nose→fwd mismatch: $(M * [0, 0, -1]) vs $fwd_scene")
    isapprox(M * [0, 1, 0], [0, 0, 1]; atol = 0.05) ||
        error("asset-up→scene-up mismatch: $(M * [0, 1, 0])")
end

# =============================== screen + sky + record ===========================================
screen = OM.Screen(scene)
OM.push_environment_image!(screen, SKY; intensity = 0.7)   # stashed; applied at author time
logmsg("first colorbuffer (stage composition + shader warmup + preroll)...")
t0 = time()
img = MK.colorbuffer(screen)
logmsg("frame 1 in $(round(time() - t0; digits = 1))s: size=$(size(img))")
for (tgt, _, _, _, _, _) in WHEELS                        # MUST_EXIST re-probe (throws if wrong)
    OM.bind_usd!(car, tgt, wobs[tgt])
end

idx_of(t) = clamp(searchsortedfirst(sol.t, t), 1, N)
still_prefix = "country_conceptcar_" * (SHOW_MBC ? "mbc_" : "") * (CAR_GLASS ? "glass_" : "") * "still_"
for (tag, i) in [("t02", idx_of(2.0)), ("t25", idx_of(25.0)), ("t45", idx_of(45.0))]
    set_frame!(i)
    OM.reset_accumulation!(screen)                        # teleport → flush reprojection history
    MK.colorbuffer(screen); MK.colorbuffer(screen)
    OM.PNGFiles.save(joinpath(DATA, still_prefix * tag * ".png"), MK.colorbuffer(screen))
    logmsg("still $tag saved")
end

set_frame!(idx_of(10.0)); OM.reset_accumulation!(screen); MK.colorbuffer(screen)
tprobe = time()
NP = 48
for i in idx_of(10.0):idx_of(10.0)+NP-1
    set_frame!(i)
    MK.colorbuffer(screen)
end
fps = NP / (time() - tprobe)
logmsg("TIMING: $(round(fps; digits = 2)) fps → full $(N)-frame video ≈ $(round(N / fps / 60; digits = 1)) min")

if get(ENV, "COUNTRY_RECORD", "1") == "1"
    set_frame!(1); OM.reset_accumulation!(screen); MK.colorbuffer(screen)   # re-converge at t=0
    H, W = size(MK.colorbuffer(screen))
    logmsg("recording $N frames at $(W)x$(H) -> $VIDEO_OUT ...")
    using FFMPEG_jll
    t_rec = time()
    FFMPEG_jll.ffmpeg() do exe
        open(pipeline(`$exe -y -f rawvideo -pixel_format rgba -video_size $(W)x$(H) -framerate 24 -i pipe:0 -c:v libx264 -preset medium -crf 20 -pix_fmt yuv420p $VIDEO_OUT`;
                      stdout = devnull, stderr = devnull), "w") do vio
            for i in 1:N
                set_frame!(i)
                write(vio, permutedims(MK.colorbuffer(screen)))   # [H,W] top-left → scanlines
                i % 240 == 0 && logmsg("frame $i/$N ($(round(i / (time() - t_rec); digits = 2)) fps)")
            end
        end
    end
    logmsg("VIDEO SAVED $VIDEO_OUT ($(round(filesize(VIDEO_OUT) / 1e6; digits = 1)) MB)")
end

close(screen)
logmsg("DONE")
