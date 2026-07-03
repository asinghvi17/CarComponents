# belgian_road_conceptcar_dashboard.jl — the belgian_road_animation.jl dashboard (3D render +
# suspension telemetry panels) with the 3D viewport RAYTRACED: the multibody primitive car is
# swapped for the NVIDIA ConceptCar01 usdplot (ghost-glass body over the live multibody
# suspension skeleton), and the LScene is embedded via OmniverseMakie.replace_scene! — GLMakie
# draws the figure + 2D axes, ovrtx RTX draws the 3D scene into the LScene's viewport rectangle.
#
# Run (GPU box):
#   flock /tmp/omniversemakie-gpu.lock env DISPLAY=:0 JULIA_CUDA_USE_COMPAT=false \
#     OVRTX_LIBRARY_PATH=<...>/ovrtx/bin/libovrtx-dynamic.so \
#     julia +dyad-3.2.0-next.87 --project=. scripts/belgian_road_conceptcar_dashboard.jl
#
# Env knobs:
#   BELGIAN_SOLVE_ONLY  "1" to stop after caching the solve (no GPU needed)
#   BELGIAN_RECORD      "0" to skip the video (preview stills + timing only)
#   CAR_GLASS           "0" for the opaque paint body (default 1 = ghost-glass shells)
#   GLASS_OPACITY       ghost-shell opacity per surface (default 0.04)
#   CONCEPTCAR_USD / SKY_HDR   asset overrides (same as the country script)
#
# The car model is the SAME exact-aspect-ratio geometry as the country-road ConceptCar video:
# axle_spacing / wheel_radius are the ConceptCar scaled to the sim's 1.48 m contact track
# (s = 0.79494), so the asset needs one uniform scale and wheelbase AND track land exactly.
#
# Load-bearing conventions (inherited from country_road_conceptcar_animation.jl — see its header):
#   • Scene is Z-UP; the MBC plots are authored in SIM (Y-up) coords and the LScene carries a
#     single Q_X90 root rotation.  The camera + lights live in world space (s3 mapping).
#   • usdplot up = :z; car rotation = rot_to_quat(R_sim)·Q_BASE (scene already carries Q_X90).
#   • Recording follows the replace_scene! docstring recipe (OmniverseMakie ≥ fefb7e4): stop the
#     GLMakie renderloop (close_after_renderloop = false!), then OM.record_frame!(session; ticks)
#     per frame — each tick is one synchronous GLMakie.colorbuffer (render_tick → sync →
#     steps_per_tick RTX steps → blit → composite), so a frame = ticks × steps_per_tick samples.

using CarComponents
using ModelingToolkit
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using ADTypes
using SciMLBase
using MultibodyComponents
using DyadCompilerPasses
using GLMakie
using Serialization
using LinearAlgebra
import OmniverseMakie as OM
const MK = OM.Makie

const ROOT = dirname(@__DIR__)
const DATA = joinpath(@__DIR__, "data")                    # gitignored
mkpath(DATA)
const CAR = get(ENV, "CONCEPTCAR_USD",
    "/home/juliahub/temp/digital-twins-for-fluid-simulation/stages/layers/CarHero/ConceptCar/ConceptCar01_Adjust.usd")
const SKY = get(ENV, "SKY_HDR",
    "/home/juliahub/.cache/packman/chk/usd.py312.manylinux_2_35_x86_64.stock.release/0.25.02.kit.8-gl.16788+c1c423f2/resources/Lights/table_mountain.hdr")
const SOL_CACHE  = joinpath(DATA, "belgian_conceptcar_sol.jls")
const CAR_GLASS  = get(ENV, "CAR_GLASS", "1") == "1"
const GLASS_OPACITY = parse(Float64, get(ENV, "GLASS_OPACITY", "0.04"))
const VIDEO_OUT = joinpath(ROOT, "assets",
    "belgian_road_conceptcar_dashboard" * (CAR_GLASS ? "_glass" : "") * "_rt.mp4")
const G0 = 9.80665
logmsg(m) = (println("[", round(time(); digits = 1), "] ", m); flush(stdout))

# ConceptCar-ratio geometry scaled to the 1.48 m contact track (same values the country-road
# model uses — ConceptCar wheelbase 2.9538 m × 0.79494, mean wheel radius × 0.79494).
const TARGET_AXLE_SPACING = 2.3481016274718134
const TARGET_WHEEL_RADIUS = 0.30760131319880757

# ===================================== model + solve =============================================
function compile_belgian_model(road_surface)
    @named model = CarComponents.ControlledBelgianRoadStraightLineCar(
        wheel_elastic_contact = true,
        axle_spacing = TARGET_AXLE_SPACING,
        wheel_radius = TARGET_WHEEL_RADIUS,
        steer_limit = 0.25,
        heading_gain = 0.8,
        lateral_gain = 0.25,
        road_surface = road_surface,
    )
    reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(
        inline_linear_sccs = true,
        analytical_linear_scc_limit = 1,
    )
    sys = ModelingToolkit.mtkcompile(model;
        additional_passes = [], reassemble_alg,
        optimize = [DyadCompilerPasses.LDIV_RULE])
    return model, sys
end

function belgian_initial_conditions(sys, road_surface; speed, baseline = nothing, start_x = 1.0)
    corners = [sys.excited_suspension_fr, sys.excited_suspension_fl,
               sys.excited_suspension_br, sys.excited_suspension_bl]
    wheel_dirs = [1.0, -1.0, 1.0, -1.0]
    defs = Pair[]
    for (corner, dir) in zip(corners, wheel_dirs)
        push!(defs, corner.wheel.wheeljoint.v_small => 1e-3)
        push!(defs, corner.suspension.ks => 5 * 44000)
        push!(defs, corner.suspension.cs => 5 * 4000)
        push!(defs, corner.suspension.r2.phi => 5.932380614359173)
        push!(defs, corner.wheel_rotation.phi => 0.0)
        push!(defs, corner.wheel_rotation.w => -dir * speed / TARGET_WHEEL_RADIUS)
    end
    # 0.193 was the settled chassis height for the r=0.2 default wheel; the bigger wheel lifts
    # the chassis by ΔR.  Only an init guess — the static settle run refines it.
    body_y0 = baseline === nothing ?
        0.193 + (TARGET_WHEEL_RADIUS - 0.2) + road_surface(start_x, 0.0) : baseline.body_y
    append!(defs, [
        sys.back_front.body.r_0[1] => start_x,
        sys.back_front.body.r_0[2] => body_y0,
        sys.back_front.body.r_0[3] => 0.0,
        sys.back_front.body.v_0[1] => speed,
        sys.back_front.body.v_0[2] => 0.0,
        sys.back_front.body.v_0[3] => 0.0,
        sys.back_front.body.phi[2] => 0.0,
        sys.excited_suspension_fr.steer_rotation.phi => 0.0,
        sys.excited_suspension_fl.steer_rotation.phi => 0.0,
        sys.excited_suspension_fr.steering_position.w => 0.0,
        sys.excited_suspension_fl.steering_position.w => 0.0,
        sys.world.render => true,
        sys.world.axis_radius => 0.01,
    ])
    if baseline !== nothing
        append!(defs, [
            sys.excited_suspension_fr.suspension.springdamper.s => baseline.fr_s,
            sys.excited_suspension_fl.suspension.springdamper.s => baseline.fl_s,
            sys.excited_suspension_br.suspension.springdamper.s => baseline.br_s,
            sys.excited_suspension_bl.suspension.springdamper.s => baseline.bl_s,
        ])
    end
    return defs
end

# Same no-init workaround as the country script: the symbolic initialization-problem construction
# hangs on this machine (SymbolicUtils hashconsing blowup on controller-feedback models), so seed
# the algebraic unknowns and let BrownFullBasicInit() initialize numerically at solve time.
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

function simulate_noinit(sys, road_surface; speed, tstop, save_dt, baseline = nothing, start_x = 1.0)
    defs = belgian_initial_conditions(sys, road_surface; speed, baseline, start_x)
    append!(defs, algebraic_seeds(sys, 0.0))
    prob = ODEProblem(sys, defs, (0.0, tstop);
        optimize = :basic, saveat = save_dt,
        warn_initialize_determined = false, build_initializeprob = false)
    sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff());
        initializealg = BrownFullBasicInit(), abstol = 1e-6, reltol = 1e-6, dtmax = 0.01)
    SciMLBase.successful_retcode(sol) || error("solve failed: $(sol.retcode), t_end=$(sol.t[end])")
    return sol
end

function solve_or_load()
    road_surface = CarComponents.belgian_block_surface_interpolator()
    if isfile(SOL_CACHE)
        logmsg("loading solve cache $SOL_CACHE ...")
        ck = deserialize(SOL_CACHE)
        return ck.sol, ck.baseline, road_surface
    end
    logmsg("compiling Belgian-road model (exact ConceptCar-ratio geometry)...")
    model, sys = compile_belgian_model(road_surface)
    logmsg("compiled: $(length(unknowns(sys))) unknowns; settling static baseline (4 s)...")
    bsol = simulate_noinit(sys, road_surface; speed = 0.0, tstop = 4.0, save_dt = 1 / 60)
    tb = bsol.t[end]
    wx = Dict(c => bsol(tb; idxs = getproperty(sys, Symbol("wheel_position_", c))) for c in (:fr, :fl, :br, :bl))
    wz = Dict(c => bsol(tb; idxs = getproperty(sys, Symbol("wheel_lateral_position_", c))) for c in (:fr, :fl, :br, :bl))
    body_x = bsol(tb; idxs = sys.back_front.body.r_0[1])
    body_z = bsol(tb; idxs = sys.back_front.body.r_0[3])
    baseline = (
        body_y = bsol(tb; idxs = sys.back_front.body.r_0[2]),
        body_ay = bsol(tb; idxs = sys.back_front.body.a_0[2]),
        body_x = body_x, body_z = body_z,
        road_body_center = road_surface(body_x, body_z),
        road_wheel_average = sum(road_surface(wx[c], wz[c]) for c in (:fr, :fl, :br, :bl)) / 4,
        fr_s = bsol(tb; idxs = sys.excited_suspension_fr.suspension.springdamper.s),
        fl_s = bsol(tb; idxs = sys.excited_suspension_fl.suspension.springdamper.s),
        br_s = bsol(tb; idxs = sys.excited_suspension_br.suspension.springdamper.s),
        bl_s = bsol(tb; idxs = sys.excited_suspension_bl.suspension.springdamper.s),
    )
    logmsg("moving run (0.8 m/s, 12 s)...")
    sol = simulate_noinit(sys, road_surface; speed = 0.8, tstop = 12.0, save_dt = 1 / 60,
        baseline = baseline)
    serialize(SOL_CACHE, (; sol, baseline))
    logmsg("solve cached: $SOL_CACHE ($(round(filesize(SOL_CACHE) / 1e6; digits = 1)) MB)")
    return sol, baseline, road_surface
end

sol, baseline, road_surface = solve_or_load()
if get(ENV, "BELGIAN_SOLVE_ONLY", "0") == "1"
    logmsg("BELGIAN_SOLVE_ONLY=1 — solve cached, exiting before render.")
    exit(0)
end
model, sys = compile_belgian_model(road_surface)
logmsg("model ready: $(length(sol.t)) samples, t=[$(sol.t[1]), $(sol.t[end])]")

isfile(CAR) || error("ConceptCar asset missing: $CAR (set CONCEPTCAR_USD)")
isfile(SKY) || error("sky HDR missing: $SKY (set SKY_HDR)")

# ============================ per-frame signals (extracted once) =================================
series(var) = collect(sol[var])
fa = (x = series(sys.back_front.frame_a.r_0[1]),           # front-axle center (anchor)
      y = series(sys.back_front.frame_a.r_0[2]),
      z = series(sys.back_front.frame_a.r_0[3]))
fb = (x = series(sys.back_front.frame_b.r_0[1]),           # rear-axle center
      y = series(sys.back_front.frame_b.r_0[2]),
      z = series(sys.back_front.frame_b.r_0[3]))
axR = (x = series(sys.front_axle.frame_a.r_0[1]),          # front axle RIGHT end (fr)
       y = series(sys.front_axle.frame_a.r_0[2]),
       z = series(sys.front_axle.frame_a.r_0[3]))
axL = (x = series(sys.front_axle.frame_b.r_0[1]),          # front axle LEFT end (fl)
       y = series(sys.front_axle.frame_b.r_0[2]),
       z = series(sys.front_axle.frame_b.r_0[3]))
steer = series(sys.controller.steer_angle)
s_fr = series(sys.excited_suspension_fr.suspension.springdamper.s)
s_fl = series(sys.excited_suspension_fl.suspension.springdamper.s)
s_br = series(sys.excited_suspension_br.suspension.springdamper.s)
s_bl = series(sys.excited_suspension_bl.suspension.springdamper.s)
N = length(sol.t)

arc = zeros(N)                                             # distance rolled, for wheel spin
for i in 2:N
    arc[i] = arc[i-1] + hypot(fa.x[i] - fa.x[i-1], fa.z[i] - fa.z[i-1])
end

# steer-sign calibration (straight-line hold → tiny angles; sign is cosmetic here)
fwdx = fa.x .- fb.x; fwdz = fa.z .- fb.z
lrate = zeros(N)
for i in 1:N-1
    lrate[i] = fwdz[i] * fwdx[i+1] - fwdx[i] * fwdz[i+1]
end
const K_STEER = sum(steer .* lrate) >= 0 ? 1.0f0 : -1.0f0

# ====================== ConceptCar01 geometry constants (all measured) ===========================
const S_SCALE = TARGET_AXLE_SPACING / 2.9538               # ≈ 0.79494
const S_TOT   = S_SCALE * 0.01                             # asset cm → world m
const R_VIS_F = 0.3827 * S_SCALE                           # rendered wheel radii [m]
const R_VIS_B = 0.3900 * S_SCALE
const C_FWD   = 0.0936 * S_SCALE                           # asset origin ahead of the front axle
δ = fa.y[1] - road_surface(fa.x[1], fa.z[1])               # settled front-axle height above road
logmsg("ride height: δ=$(round(δ; digits = 4)) m")

function body_pose(i)                                      # sim-frame R (columns fwd/up/right), O
    X = normalize([fwdx[i], fa.y[i] - fb.y[i], fwdz[i]])
    lat = [axR.x[i] - axL.x[i], axR.y[i] - axL.y[i], axR.z[i] - axL.z[i]]
    Z = normalize(lat .- (dot(lat, X) .* X))
    Y = cross(Z, X)
    R = [X Y Z]
    O = [fa.x[i], fa.y[i], fa.z[i]] .+ R * [C_FWD, -δ, 0.0]
    return R, O
end

function rot_to_quat(R)                                    # self-checked vs Makie's convention
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
function _noise_octave(nextf, px, cells)                   # periodic bilinear value noise
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
    lo = _noise_octave(nextf, px, 4)                       # broad dry patches
    img = map(n, lo) do t, d
        g1 = (0.16f0, 0.25f0, 0.09f0); g2 = (0.34f0, 0.44f0, 0.16f0); dry = (0.38f0, 0.34f0, 0.18f0)
        r = g1[1] + (g2[1] - g1[1]) * t; g = g1[2] + (g2[2] - g1[2]) * t; b = g1[3] + (g2[3] - g1[3]) * t
        w = clamp((d - 0.62f0) * 2.2f0, 0f0, 0.28f0)
        MK.RGBf(r + (dry[1] - r) * w, g + (dry[2] - g) * w, b + (dry[3] - b) * w)
    end
    OM.PNGFiles.save(path, img)
end
isfile(GRASS_PNG) || make_grass!(GRASS_PNG)

# Cobblestone albedo straight from the CRG height field: seams (local lows) darken, crowns
# lighten — the texture is spatially registered with the geometry the suspension actually feels.
const COBBLE_PNG = joinpath(DATA, "belgian_cobble_albedo.png")
function make_cobble!(path, road_surface; x0, x1, z0, z1, ppm = 200)
    nx = clamp(round(Int, (x1 - x0) * ppm), 512, 2800)
    nz = clamp(round(Int, (z1 - z0) * ppm), 128, 1024)
    Hf = [Float32(road_surface(x0 + (i - 0.5) * (x1 - x0) / nx, z0 + (j - 0.5) * (z1 - z0) / nz))
          for i in 1:nx, j in 1:nz]
    finite = filter(isfinite, vec(Hf))                     # the CRG grid has NaN dropout cells
    lo, hi = extrema(finite)
    rng = max(hi - lo, 1.0f-6)
    nextf = _lcg(UInt64(0x4d595df4d0f33173))
    spek = _noise_octave(nextf, 1024, 256)                 # fine granite speckle (tiled sampling)
    img = Matrix{MK.RGBf}(undef, nx, nz)
    for j in 1:nz, i in 1:nx
        t = isfinite(Hf[i, j]) ? (Hf[i, j] - lo) / rng : 0.35f0   # NaN cells → mid-tone
        s = spek[(i - 1) % 1024 + 1, (j - 1) % 1024 + 1]
        g = clamp(0.15f0 + 0.36f0 * t + 0.10f0 * (s - 0.5f0), 0.05f0, 0.8f0)
        img[i, j] = MK.RGBf(g * 1.03f0, g * 0.99f0, g * 0.94f0)   # warm-gray granite
    end
    OM.PNGFiles.save(path, img)
    return nothing
end

# road extents from the CRG data itself (the interpolator extrapolates flat beyond them)
xaxis_crg, zaxis_crg, _H = CarComponents.opencrg_surface_data(CarComponents.belgian_block_crg_path())
const RX0, RX1 = max(-0.5, first(xaxis_crg)), min(12.0, last(xaxis_crg))
const RZ0, RZ1 = first(zaxis_crg), last(zaxis_crg)
logmsg("belgian CRG extent: x ∈ [$(round(RX0; digits=2)), $(round(RX1; digits=2))], z ∈ [$(round(RZ0; digits=2)), $(round(RZ1; digits=2))]")
isfile(COBBLE_PNG) || make_cobble!(COBBLE_PNG, road_surface; x0 = RX0, x1 = RX1, z0 = RZ0, z1 = RZ1)

# =========================== dashboard data (2D telemetry series) ================================
function dashboard_data(sys, sol, road_surface, baseline)
    body_x = sol[sys.back_front.body.r_0[1]]
    body_z = sol[sys.back_front.body.r_0[3]]
    fr_x = sol[sys.wheel_position_fr]; fl_x = sol[sys.wheel_position_fl]
    br_x = sol[sys.wheel_position_br]; bl_x = sol[sys.wheel_position_bl]
    fr_z = sol[sys.wheel_lateral_position_fr]; fl_z = sol[sys.wheel_lateral_position_fl]
    br_z = sol[sys.wheel_lateral_position_br]; bl_z = sol[sys.wheel_lateral_position_bl]
    road_body_center_mm = [1000 * (road_surface(body_x[i], body_z[i]) - baseline.road_body_center)
                           for i in eachindex(body_x)]
    road_wheel_average_mm = [1000 * ((road_surface(fr_x[i], fr_z[i]) + road_surface(fl_x[i], fl_z[i]) +
                                      road_surface(br_x[i], br_z[i]) + road_surface(bl_x[i], bl_z[i])) / 4 -
                                     baseline.road_wheel_average) for i in eachindex(body_x)]
    return (
        body_x = body_x, body_z = body_z,
        body_ay_g = (sol[sys.back_front.body.a_0[2]] .- baseline.body_ay) ./ G0,
        heading = -sol[sys.back_front.body.phi[2]],
        steer = sol[sys.controller.steer_angle],
        road_body_center_mm = road_body_center_mm,
        road_wheel_average_mm = road_wheel_average_mm,
        fr_s = 1000 .* (sol[sys.excited_suspension_fr.suspension.springdamper.s] .- baseline.fr_s),
        fl_s = 1000 .* (sol[sys.excited_suspension_fl.suspension.springdamper.s] .- baseline.fl_s),
        br_s = 1000 .* (sol[sys.excited_suspension_br.suspension.springdamper.s] .- baseline.br_s),
        bl_s = 1000 .* (sol[sys.excited_suspension_bl.suspension.springdamper.s] .- baseline.bl_s),
    )
end

function padded_limits(values; pad_fraction = 0.08, min_span = 1e-6)
    lo, hi = extrema(values)
    span = max(hi - lo, min_span)
    return lo - pad_fraction * span, hi + pad_fraction * span
end

# ============================== figure + 3D scene (Z-up) =========================================
OM.activate!(warmup = 24, accumulate_across_frames = true, background = :domelight)
GLMakie.activate!()

lights = MK.AbstractLight[
    MK.DirectionalLight(MK.RGBf(1.8, 1.75, 1.6), s3(-0.45, -1.0, -0.3)),   # sun (dome adds ambient)
]

fig, time_obs, ls = MultibodyComponents.render(model, sol, sol.t[1];
    slider = false, show_axis = false, size = (1920, 1080), lights = lights)
scene = ls.scene
MK.rotate!(scene, Q_X90)                                   # SIM (Y-up) content → Z-up world
axplots = filter(p -> p isa MK.Axis3D, copy(scene.plots))
foreach(p -> delete!(scene, p), axplots)
# drop the red plots (SlippingWheel cylinders — the ConceptCar wheels replace them — and the
# world x-axis arrow); force the kept hardware bright + opaque (dyad gray is 30% alpha).
redplots = filter(copy(scene.plots)) do p
    c = try p.color[] catch; nothing end
    c isa MK.Colors.Colorant && Float32(MK.Colors.red(c)) > 0.9f0 &&
        Float32(MK.Colors.green(c)) < 0.15f0 && Float32(MK.Colors.blue(c)) < 0.15f0
end
foreach(p -> delete!(scene, p), redplots)
recolored = Ref(0)
for p in copy(scene.plots)
    c = try p.color[] catch; nothing end
    c isa MK.Colors.Colorant || continue
    r = Float32(MK.Colors.red(c)); g = Float32(MK.Colors.green(c)); b = Float32(MK.Colors.blue(c))
    p.color[] = MK.RGBAf(min(1, 0.3f0 + r), min(1, 0.3f0 + g), min(1, 0.3f0 + b), 1.0f0)
    recolored[] += 1
end
logmsg("MBC scene: removed $(length(axplots)) axis + $(length(redplots)) red plots; " *
       "$(length(scene.plots)) kept ($(recolored[]) recolored opaque)")
W3(x, y, z) = MK.Point3f(x, y, z)                          # content authored in SIM coords

# static camera, the original dashboard's framing pulled back a touch (car drives through the
# shot) and tilted up for a sliver of horizon.  Applied here AND per frame in set_frame! — the
# MBC render / Camera3D re-fit the camera around display time, so a single pre-display
# update_cam! silently loses (probe-verified).
set_camera!() = MK.update_cam!(scene, s3(5.5, 3.4, 7.0), s3(5.5, 0.75, 0.0), MK.Vec3f(0, 0, 1))
set_camera!()

# cobbled road strip with the CRG-derived texture (one UV tile over the whole strip)
function road_mesh_uv(; nx = 560, nz = 64)
    xs = range(RX0, RX1; length = nx)
    zs = range(RZ0, RZ1; length = nz)
    points = MK.Point3f[]; uv = MK.Vec2f[]
    for j in eachindex(zs), i in eachindex(xs)
        h = road_surface(xs[i], zs[j])
        isfinite(h) || (h = 0.0)                           # NaN dropout cells in the CRG grid
        push!(points, W3(xs[i], h, zs[j]))
        push!(uv, MK.Vec2f((xs[i] - RX0) / (RX1 - RX0), (zs[j] - RZ0) / (RZ1 - RZ0)))
    end
    faces = GLMakie.GeometryBasics.GLTriangleFace[]
    for j in 1:(nz - 1), i in 1:(nx - 1)
        a = (j - 1) * nx + i; b = a + 1; c = j * nx + i; d = c + 1
        push!(faces, GLMakie.GeometryBasics.GLTriangleFace(a, c, b))
        push!(faces, GLMakie.GeometryBasics.GLTriangleFace(b, c, d))
    end
    return GLMakie.GeometryBasics.Mesh(points, faces; uv = uv)
end
MK.mesh!(scene, road_mesh_uv();
         color = :white, material = (; base_color_texture = COBBLE_PNG,
                                       roughness = 0.85f0, metallic = 0.0f0))

# grass ground: one textured quad just under the road, UVs tiled every ~5 m.  The plane must sit
# below the DEEPEST cobble hollow or it pokes through the road as green patches.
hmin_road = minimum(filter(isfinite,
    [road_surface(x, z) for x in range(RX0, RX1; length = 160), z in range(RZ0, RZ1; length = 40)]))
hg = hmin_road - 0.03
gx0, gx1 = RX0 - 60, RX1 + 60
gz0, gz1 = RZ0 - 60, RZ1 + 60
ntx = round((gx1 - gx0) / 5); nty = round((gz1 - gz0) / 5)
gpts = [W3(gx0, hg, gz0), W3(gx1, hg, gz0), W3(gx1, hg, gz1), W3(gx0, hg, gz1)]
guv  = MK.Vec2f[(0, 0), (ntx, 0), (ntx, nty), (0, nty)]
gnrm = [MK.Vec3f(0, 1, 0) for _ in 1:4]                    # sim-frame up
gfcs = GLMakie.GeometryBasics.GLTriangleFace[(1, 2, 3), (1, 3, 4)]
MK.mesh!(scene, GLMakie.GeometryBasics.Mesh(gpts, gfcs; uv = guv, normal = gnrm); color = :white,
         material = (; base_color_texture = GRASS_PNG, roughness = 0.95f0, metallic = 0.0f0))

# center reference line + driven front-axle path (thin world-space tubes, meters)
road_line_x = collect(range(max(RX0, 0.0), RX1; length = 200))
MK.lines!(scene, [W3(x, road_surface(x, 0.0) + 0.025, 0.0) for x in road_line_x];
          color = :deepskyblue, linewidth = 0.04)
MK.lines!(scene, [W3(fa.x[i], road_surface(fa.x[i], fa.z[i]) + 0.04, fa.z[i]) for i in 1:N];
          color = :orange, linewidth = 0.035)

# ================================== the car ======================================================
CAR_EFFECTIVE = CAR
if CAR_GLASS
    wrap = joinpath(DATA, "conceptcar_glass_wrapper_belgian.usda")
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
    logmsg("glass wrapper written: $wrap (opacity $(GLASS_OPACITY))")
end
car = OM.usdplot!(scene, CAR_EFFECTIVE; up = :z,
    bbox = MK.Rect3f(MK.Point3f(-111.4, -0.1, -100.7), MK.Vec3f(222.7, 137.4, 505.0)))
MK.scale!(car, S_TOT, S_TOT, S_TOT)

# Wheel groups /root/Wheels/wheel_{DF,PF,DB,PB} (D=driver=asset−x=fl/bl, P=passenger=fr/br);
# parent frame is cm-scale, x=LEFT, y=up, z=nose: +R_x = forward roll, +R_y = steer left.
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

# ========================== 2D telemetry panels (right column) ===================================
data = dashboard_data(sys, sol, road_surface, baseline)
xplot_hi = ceil(maximum(data.body_x) * 5) / 5 + 0.2

side = fig[1, 2] = GLMakie.GridLayout()
GLMakie.colsize!(fig.layout, 1, GLMakie.Relative(0.62))
GLMakie.colsize!(fig.layout, 2, GLMakie.Relative(0.38))
GLMakie.Label(fig[0, :],
    "Controlled straight-line car on Belgian-block road: raytraced render + suspension plots";
    fontsize = 26, tellwidth = false)

status = GLMakie.Observable("t = 0.00 s    x = 0.00 m    ay = 0.00 g")
GLMakie.Label(side[1, 1], status; fontsize = 18, tellwidth = false, halign = :left)

ax1 = GLMakie.Axis(side[2, 1]; xlabel = "distance along road x [m]",
    ylabel = "relative road height [mm]", title = "2D road input relative to standing still")
GLMakie.lines!(ax1, data.body_x, data.road_wheel_average_mm; color = :black, linewidth = 2,
    label = "average road under wheels")
GLMakie.lines!(ax1, data.body_x, data.road_body_center_mm; color = :gray45, linestyle = :dash,
    linewidth = 2, label = "road at chassis X-Z")
GLMakie.axislegend(ax1; position = :lb, labelsize = 12)
GLMakie.xlims!(ax1, 0, xplot_hi)
GLMakie.ylims!(ax1, padded_limits(vcat(data.road_wheel_average_mm, data.road_body_center_mm))...)

ax2 = GLMakie.Axis(side[3, 1]; xlabel = "distance along road x [m]",
    ylabel = "vertical acceleration [g]", title = "Chassis center vertical acceleration")
GLMakie.hlines!(ax2, [0.0]; color = :gray70, linewidth = 1)
GLMakie.lines!(ax2, data.body_x, data.body_ay_g; color = :royalblue, linewidth = 2, label = "chassis ay")
GLMakie.axislegend(ax2; position = :lb, labelsize = 12)
GLMakie.xlims!(ax2, 0, xplot_hi)
GLMakie.ylims!(ax2, padded_limits(data.body_ay_g)...)

ax3 = GLMakie.Axis(side[4, 1]; xlabel = "time [s]",
    ylabel = "spring travel relative to static [mm]", title = "Suspension travel check")
GLMakie.lines!(ax3, sol.t, data.fr_s; color = :red, linewidth = 2, label = "FR")
GLMakie.lines!(ax3, sol.t, data.fl_s; color = :magenta, linewidth = 2, label = "FL")
GLMakie.lines!(ax3, sol.t, data.br_s; color = :green, linewidth = 2, label = "BR")
GLMakie.lines!(ax3, sol.t, data.bl_s; color = :cyan4, linewidth = 2, label = "BL")
GLMakie.axislegend(ax3; position = :lb, nbanks = 2, labelsize = 12)
GLMakie.xlims!(ax3, sol.t[1], sol.t[end])
# explicit ylims on ax3/ax4: the t_marker vlines! poison Makie's y-autolimits (axis falls back
# to the (0,10) default and the data clips out of view)
GLMakie.ylims!(ax3, padded_limits(vcat(data.fr_s, data.fl_s, data.br_s, data.bl_s))...)

ax4 = GLMakie.Axis(side[5, 1]; xlabel = "time [s]",
    ylabel = "z [m] / angles [rad]", title = "Straight-line tracking and steering effort")
GLMakie.lines!(ax4, sol.t, data.body_z; color = :orange, linewidth = 2, label = "lateral z")
GLMakie.lines!(ax4, sol.t, data.heading; color = :purple, linewidth = 2, label = "heading")
GLMakie.lines!(ax4, sol.t, data.steer; color = :brown, linewidth = 2, label = "steer")
GLMakie.axislegend(ax4; position = :lb, nbanks = 3, labelsize = 12)
GLMakie.xlims!(ax4, sol.t[1], sol.t[end])
GLMakie.ylims!(ax4, padded_limits(vcat(data.body_z, data.heading, data.steer))...)

x_marker = GLMakie.Observable([data.body_x[1]])
t_marker = GLMakie.Observable([sol.t[1]])
road_marker = GLMakie.Observable([data.road_wheel_average_mm[1]])
accel_marker = GLMakie.Observable([data.body_ay_g[1]])
GLMakie.vlines!(ax1, x_marker; color = :gray30, linewidth = 2)
GLMakie.scatter!(ax1, x_marker, road_marker; color = :black, markersize = 12)
GLMakie.vlines!(ax2, x_marker; color = :gray30, linewidth = 2)
GLMakie.scatter!(ax2, x_marker, accel_marker; color = :royalblue, markersize = 12)
GLMakie.vlines!(ax3, t_marker; color = :gray30, linewidth = 2)
GLMakie.vlines!(ax4, t_marker; color = :gray30, linewidth = 2)

function update_markers!(time)
    x = sol(time; idxs = sys.back_front.body.r_0[1])
    z = sol(time; idxs = sys.back_front.body.r_0[3])
    wfr_x = sol(time; idxs = sys.wheel_position_fr); wfl_x = sol(time; idxs = sys.wheel_position_fl)
    wbr_x = sol(time; idxs = sys.wheel_position_br); wbl_x = sol(time; idxs = sys.wheel_position_bl)
    wfr_z = sol(time; idxs = sys.wheel_lateral_position_fr); wfl_z = sol(time; idxs = sys.wheel_lateral_position_fl)
    wbr_z = sol(time; idxs = sys.wheel_lateral_position_br); wbl_z = sol(time; idxs = sys.wheel_lateral_position_bl)
    road_average_mm = 1000 * ((road_surface(wfr_x, wfr_z) + road_surface(wfl_x, wfl_z) +
                               road_surface(wbr_x, wbr_z) + road_surface(wbl_x, wbl_z)) / 4 -
                              baseline.road_wheel_average)
    ay_g = (sol(time; idxs = sys.back_front.body.a_0[2]) - baseline.body_ay) / G0
    x_marker[] = [x]; t_marker[] = [time]
    road_marker[] = [road_average_mm]; accel_marker[] = [ay_g]
    status[] = "t = $(round(time, digits = 2)) s    x = $(round(x, digits = 2)) m    " *
               "z = $(round(z, digits = 3)) m    ay = $(round(ay_g, digits = 2)) g"
    return nothing
end

# =============================== per-frame update ================================================
function set_frame!(i)
    R, O = body_pose(i)
    MK.rotate!(car, rot_to_quat(R) * Q_BASE)               # scene root already carries Q_X90
    MK.translate!(car, MK.Vec3f(W3(O...)))
    time_obs[] = sol.t[i]                                  # drive the MBC plots
    update_markers!(sol.t[i])
    for (tgt, hub, steerable, L, rvis, sser) in WHEELS
        dy = Float32(-(sser[i] - sser[1]) / S_TOT)         # spring shortens → wheel UP (parent cm)
        st = steerable ? K_STEER * Float32(steer[i]) : 0.0f0
        sp = Float32(mod2pi(arc[i] / rvis))                # +R_x = forward roll
        wobs[tgt][] = wheel_mat(hub, dy, st, sp, L)
    end
    set_camera!()                                          # win the last-write race (see above)
    return nothing
end
set_frame!(1)

# fail-fast: composed rotation maps asset nose(−z)→scene fwd and asset up(y)→scene +z
let (R1, _) = body_pose(1)
    M = Matrix(MK.rotationmatrix4(Q_X90 * rot_to_quat(R1) * Q_BASE))[1:3, 1:3]
    h1 = atan(fwdz[1], fwdx[1])
    fwd_scene = [cos(h1), -sin(h1), 0.0]
    isapprox(M * [0, 0, -1], fwd_scene; atol = 0.05) ||
        error("nose→fwd mismatch: $(M * [0, 0, -1]) vs $fwd_scene")
    isapprox(M * [0, 1, 0], [0, 0, 1]; atol = 0.05) ||
        error("asset-up→scene-up mismatch: $(M * [0, 1, 0])")
end

# ====================== display + replace_scene! (hybrid) + record ===============================
glscr = GLMakie.Screen(; visible = false, px_per_unit = 1, scalefactor = 1)
display(glscr, fig.scene)
GLMakie.colorbuffer(glscr)                                 # layout pass (viewport rects settle)
# The replace_scene! recording recipe: STOP the async renderloop (the embedded blit re-dirties
# the scene every tick, so a running loop self-sustains and starves libuv — pipe writes hang) and
# drive frames synchronously via record_frame!.  close_after_renderloop must be false (the
# default CLOSES the screen).
GLMakie.stop_renderloop!(glscr; close_after_renderloop = false)
logmsg("LScene viewport: $(ls.scene.viewport[]) (renderloop stopped, screen open=$(isopen(glscr)))")

session = OM.replace_scene!(ls; steps_per_tick = 8)
logmsg("embedded session up: $(length(session.screen.plot2robj)) plots authored, " *
       "fb_size=$(session.screen.fb_size)")
OM.push_environment_image!(session.screen, SKY; intensity = 0.7)   # live swap (post-author)

# One recorded frame = 3 synchronous host ticks (3 × steps_per_tick = 24 RTX samples: the first
# tick syncs the moved car and resets accumulation, the rest refine) — record_frame! is the
# scripted-recording companion to replace_scene!; see its docstring's recipe.
function frame!(i)
    set_frame!(i)
    return OM.record_frame!(session; ticks = 3)
end

for (tgt, _, _, _, _, _) in WHEELS                         # MUST_EXIST re-probe (throws if wrong)
    OM.bind_usd!(car, tgt, wobs[tgt])
end

idx_of(t) = clamp(searchsortedfirst(sol.t, t), 1, N)
prefix = "belgian_dashboard_" * (CAR_GLASS ? "glass_" : "")
for (tag, i) in [("t00", 1), ("t04", idx_of(4.0)), ("t09", idx_of(9.0))]
    img = frame!(i)
    OM.PNGFiles.save(joinpath(DATA, prefix * "still_" * tag * ".png"), img)
    logmsg("still $tag saved ($(size(img)))")
end

tprobe = time()
NP = 16
for i in idx_of(6.0):idx_of(6.0)+NP-1
    frame!(i)
end
fps = NP / (time() - tprobe)

const FRAMERATE = 30
times = collect(range(sol.t[1], sol.t[end]; step = 1 / FRAMERATE))
logmsg("TIMING: $(round(fps; digits = 2)) fps → $(length(times))-frame video ≈ " *
       "$(round(length(times) / fps / 60; digits = 1)) min")

if get(ENV, "BELGIAN_RECORD", "1") == "1"
    img1 = frame!(1)
    H, W = size(img1)
    fmt = sizeof(eltype(img1)) == 4 ? "rgba" : "rgb24"     # RGBA{N0f8} vs RGB{N0f8}
    logmsg("recording $(length(times)) frames at $(W)x$(H) ($fmt) -> $VIDEO_OUT ...")
    using FFMPEG_jll
    t_rec = time()
    FFMPEG_jll.ffmpeg() do exe
        open(pipeline(`$exe -y -f rawvideo -pixel_format $fmt -video_size $(W)x$(H) -framerate $FRAMERATE -i pipe:0 -c:v libx264 -preset medium -crf 20 -pix_fmt yuv420p $VIDEO_OUT`;
                      stdout = devnull, stderr = devnull), "w") do vio
            for (k, t) in enumerate(times)
                write(vio, permutedims(frame!(idx_of(t))))   # [H,W] top-left → scanlines
                k % 90 == 0 && logmsg("frame $k/$(length(times)) ($(round(k / (time() - t_rec); digits = 2)) fps)")
            end
        end
    end
    logmsg("VIDEO SAVED $VIDEO_OUT ($(round(filesize(VIDEO_OUT) / 1e6; digits = 1)) MB)")
end

close(session)
close(glscr)
logmsg("DONE")
