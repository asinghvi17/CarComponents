# animation3.jl — raytraced ControlledFlatRoadCar video with a REAL CAR MESH
# (NVIDIA AmericanoI7EV, textured), rendered via OmniverseMakie.
#
# Same simulation as animation2.jl; the visual swaps the primitive cylinder car
# for the AmericanoI7EV OBJ meshes (body + 4 independently spinning/steering
# wheels + glass), textured with the asset's own albedo maps, and records with
# the fast "no accumulation reset" RT2 path (~4 min for 45 s @ 30 fps @ 720p).
#
# Mesh/texture assets: bug-reports/car-mesh-sourcing/AmericanoI7EV_obj_tex/
# (extracted from the DSX content pack — NVIDIA proprietary license: fine for
# internal dev/eval, do NOT redistribute; see REPORT.md there).
using Revise
using CarComponents
using ModelingToolkit
using OrdinaryDiffEqDefault
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using ADTypes
using SciMLBase
using LinearAlgebra
using MultibodyComponents
using DyadCompilerPasses
using OmniverseMakie

const OV = OmniverseMakie
const GB = OmniverseMakie.GeometryBasics
const Makie = OmniverseMakie.Makie                # transform! is not re-exported by OV
const FileIO = OmniverseMakie.Makie.FileIO        # MeshIO loader for UV'd OBJs
const Rotations = MultibodyComponents.Rotations   # not a direct project dep; reach via MBC

const ROOT    = dirname(@__DIR__)
const MESHDIR = joinpath(ROOT, "bug-reports", "car-mesh-sourcing", "AmericanoI7EV_obj_tex")
logmsg(m) = (println("[", round(time(); digits=1), "] ", m); flush(stdout))

# ======================= simulation (same as animation2.jl) ==================
@named controlled_car = ControlledFlatRoadCar(path_center_x=0.0, path_center_z=-6.0, path_radius=6.0, steer_limit=0.6, path_direction=-1.0, heading_gain=1.0, radial_gain=0.15)
reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1)
println("Compiling controlled car...")
controlled_sys = ModelingToolkit.mtkcompile(controlled_car; additional_passes = [], reassemble_alg, optimize = [DyadCompilerPasses.LDIV_RULE])
println("compiled unknowns=", length(unknowns(controlled_sys)), " equations=", length(equations(controlled_sys)))

function simulate_controlled_car(; speed=0.8, tstop=45.0, save_dt=1/30, init_steer=atan(1/6))
    C = [controlled_sys.excited_suspension_fr, controlled_sys.excited_suspension_fl,
         controlled_sys.excited_suspension_br, controlled_sys.excited_suspension_bl]
    dirs = [1.0, -1.0, 1.0, -1.0]
    defs = Pair[]
    for (c, dir) in zip(C, dirs)
        push!(defs, c.wheel.wheeljoint.v_small => 1e-3)
        push!(defs, c.suspension.ks => 5 * 44000)
        push!(defs, c.suspension.cs => 5 * 4000)
        push!(defs, c.suspension.r2.phi => 5.932380614359173)
        push!(defs, c.wheel_rotation.phi => 0.0)
        push!(defs, c.wheel_rotation.w => -dir * speed / 0.2)
        # Seed the algebraic unknowns the symbolic initialization problem would
        # otherwise solve for (we build with build_initializeprob=false below);
        # BrownFullBasicInit refines them to consistency at solve time.
        push!(defs, c.wheel.wheeljoint.delta_0[1] => 0.0)
        push!(defs, c.wheel.wheeljoint.delta_0[3] => 0.0)
        push!(defs, c.rotational_losses.w_rel => 0.0)
        # No drive torque in the model — the car coasts; tire-slip scrub (not this
        # damper) parks it by ~t=35-40s. Reduced losses kept for the tail.
        push!(defs, c.rotational_losses.d => 0.01)
    end
    append!(defs, [
        controlled_sys.back_front.body.r_0[1] => 0.0,
        controlled_sys.back_front.body.r_0[2] => 0.193,
        controlled_sys.back_front.body.r_0[3] => 0.0,
        controlled_sys.back_front.body.v_0[1] => speed,
        controlled_sys.back_front.body.v_0[2] => 0.0,
        controlled_sys.back_front.body.v_0[3] => 0.0,
        controlled_sys.back_front.body.phi[2] => 0.0,
        controlled_sys.excited_suspension_fr.steer_rotation.phi => init_steer,
        controlled_sys.excited_suspension_fl.steer_rotation.phi => init_steer,
        controlled_sys.excited_suspension_fr.steering_position.w => 0.0,
        controlled_sys.excited_suspension_fl.steering_position.w => 0.0,
        controlled_sys.world.render => true,
        controlled_sys.world.axis_radius => 0.01,
    ])
    # NOTE: build_initializeprob=false is REQUIRED on this machine — the default
    # symbolic initialization-problem build hangs ODEProblem construction for hours
    # (SymbolicUtils hashconsing blowup; see bug-reports/init-hashconsing-hang/).
    prob = ODEProblem(controlled_sys, defs, (0.0, tstop); optimize = :basic, saveat = save_dt, build_initializeprob = false)
    solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); initializealg = BrownFullBasicInit(), abstol=1e-6, reltol=1e-6, dtmax=0.02)
end

println("Simulating controlled path follower...")
sol = simulate_controlled_car()
println("retcode=", sol.retcode, " samples=", length(sol.t), " t_end=", sol.t[end])

center_x = 0.0; center_z = -6.0; radius = 6.0

# ======================= mesh loading =========================================
# Vertex pre-transform baked ONCE at load: p_scene = s .* ( -(x+0.096), z, y ),
# s = 1/3 — rotates the asset's Z-up to the scene's Y-up, flips the nose to +X,
# centers the wheelbase at x=0 and scales wheelbase 3.0 m -> 1.0 m (= sim).
# UV'd OBJs carry face-varying uv/normal FaceViews: expand to per-vertex arrays
# (OmniverseMakie authors the `st` primvar only when length(uv)==length(points)).
const S = 1/3
@inline meshxf(p) = OV.Point3f(-S*(p[1] + 0.096), S*p[3], S*p[2])

function load_body(path)
    me  = GB.expand_faceviews(FileIO.load(path))
    fcs = GB.faces(me)
    newpos = meshxf.(GB.coordinates(me))
    uv  = GB.vertex_attributes(me).uv
    nrm = GB.normals(newpos, fcs)                # recompute smooth normals
    GB.Mesh(newpos, fcs; uv = uv, normal = nrm)
end

function load_glass(path)                        # untextured: uv unneeded
    me = GB.expand_faceviews(FileIO.load(path))
    GB.normal_mesh(meshxf.(GB.coordinates(me)), GB.faces(me))
end

# Wheels are baked at world positions in the OBJ: recenter each on its own
# transformed centroid (over the UNIQUE pre-expand verts) so it spins about its
# origin; c_i is the wheel-center offset from the body origin.
function load_wheel(path)
    m   = FileIO.load(path)
    me  = GB.expand_faceviews(m)
    fcs = GB.faces(me)
    cu  = meshxf.(GB.coordinates(m))
    c   = sum(cu) / length(cu)
    newpos = [meshxf(p) - c for p in GB.coordinates(me)]
    uv  = GB.vertex_attributes(me).uv
    nrm = GB.normals(newpos, fcs)
    GB.Mesh(newpos, fcs; uv = uv, normal = nrm), Float64[c[1], c[2], c[3]]
end

logmsg("loading UV'd OBJ meshes from $MESHDIR ...")
body_mesh  = load_body(joinpath(MESHDIR, "Americanol7EV_A1_1.obj"))
glass_mesh = load_glass(joinpath(MESHDIR, "Americanol7EV_Glass_A1_1.obj"))
wheelFR, cFR = load_wheel(joinpath(MESHDIR, "Americanol7EV_FrontWheelR_A1_1.obj"))
wheelFL, cFL = load_wheel(joinpath(MESHDIR, "Americanol7EV_FrontWheelL_A1_1.obj"))
wheelBR, cBR = load_wheel(joinpath(MESHDIR, "Americanol7EV_RearWheelR_A1_1.obj"))
wheelBL, cBL = load_wheel(joinpath(MESHDIR, "Americanol7EV_RearWheelL_A1_1.obj"))
logmsg("centroids  FR=$(round.(cFR;digits=3))  FL=$(round.(cFL;digits=3))  BR=$(round.(cBR;digits=3))  BL=$(round.(cBL;digits=3))")

# ======================= scene ================================================
# NOTE: `warmup` doubles as RTX steps/frame because OV.reset! is no-op'd for the
# recording below (accumulation carries across frames, interactive-viewport style).
OV.activate!(warmup = 4)

amb = fill(OV.RGBf(0.30, 0.33, 0.38), 1, 1)
lights = OV.AbstractLight[
    OV.EnvironmentLight(0.22, amb'),
    OV.PointLight(OV.RGBf(240, 205, 165), OV.Vec3f(center_x + 5, 16, center_z + 4)),   # warm key
    OV.PointLight(OV.RGBf(55, 70, 100),   OV.Vec3f(center_x - 9, 8,  center_z - 12)),  # cool fill
    OV.PointLight(OV.RGBf(95, 85, 105),   OV.Vec3f(center_x - 2, 10, center_z + 12)),  # subtle rim
]
fig, t, scene = MultibodyComponents.render(controlled_car, sol, sol.t[1];
    slider = false, x = 4.2, y = 6.2, z = 4.2,
    lookat = OV.Vec3f(0.0, 0.0, -6.0), up = OV.Vec3f(0.0, 1.0, 0.0),
    traces = [controlled_sys.back_front.frame_cm], size = (1280, 720), lights = lights)

# DELETE every plot render() created EXCEPT the frame_cm trace (the only Lines
# with 500 points) so only the mesh car is visible. plot.visible[]=false does
# NOT work here: OmniverseMakie authors the USD stage from scene.plots at the
# first render and initial visibility is ignored (only post-authoring changes
# propagate) — deletion pre-authoring keeps the primitives out of the stage.
npoints(p) = try length(p[1][]) catch; -1 end
function delete_primitives!(scn0)
    ndel = 0
    stack = Any[scn0]
    while !isempty(stack)
        s = pop!(stack)
        for p in collect(s.plots)                # snapshot: we mutate s.plots
            (p isa OV.Lines && npoints(p) == 500) && continue   # keep the trace
            try; delete!(s, p); ndel += 1; catch e; logmsg("  WARN delete failed on $(typeof(p)): $e"); end
        end
        append!(stack, s.children)
    end
    ndel
end
scn0 = hasproperty(scene, :scene) ? getproperty(scene, :scene) : scene
logmsg("deleted $(delete_primitives!(scn0)) primitive plots")

# ground quad + green reference circle (the controller's target path)
ground_half = radius + 8.0
ground = OV.Rect3f(OV.Vec3f(center_x - ground_half, -0.05, center_z - ground_half), OV.Vec3f(2ground_half, 0.05, 2ground_half))
OV.mesh!(scene, ground; color = :gainsboro)
let θs = range(0, 2π, length = 260)
    ref_pts = [OV.Point3f(center_x + radius*cos(a), 0.03, center_z + radius*sin(a)) for a in θs]
    OV.lines!(scene, ref_pts; color = OV.RGBf(0.12, 0.8, 0.2), linewidth = 1.2)
end

# Textured car plots. base_color_texture references the on-disk PNGs directly —
# do NOT use `color = img::Matrix` here: that path rewrites a per-plot temp PNG
# which races ovrtx's async texture reads across frames (Corrupt PNG; textures
# then survive only the first frame). `color = :white` is harmless (the bound
# diffuse_texture overrides the color constant).
bodymat  = (; base_color_texture = joinpath(MESHDIR, "body_albedo_atlas.png"), metallic = 0.1f0, roughness = 0.45f0)
wheelmat = (; base_color_texture = joinpath(MESHDIR, "wheel_albedo.png"),      metallic = 0.1f0, roughness = 0.6f0)
glassmat = (; glass = true, ior = 1.45f0, thin_walled = true)

body_plot  = OV.mesh!(scene, body_mesh;  color = :white,     material = bodymat)
glass_plot = OV.mesh!(scene, glass_mesh; color = :lightblue, material = glassmat)
fr_plot = OV.mesh!(scene, wheelFR; color = :white, material = wheelmat)
fl_plot = OV.mesh!(scene, wheelFL; color = :white, material = wheelmat)
br_plot = OV.mesh!(scene, wheelBR; color = :white, material = wheelmat)
bl_plot = OV.mesh!(scene, wheelBL; color = :white, material = wheelmat)

# ======================= per-frame pose hook ==================================
frame_cm = controlled_sys.back_front.frame_cm
h0 = get_trans(sol, frame_cm, sol.t[1])[2]       # CM height at t0 (~0.193)
o  = Float64[0.0, -h0, 0.0]                       # local offset: skin ground -> y≈0

Qf(M) = (q = Rotations.QuatRotation(M).q; OV.Quaternionf(q.v1, q.v2, q.v3, q.s))

# (plot, wheel-center offset, spin sign, suspension subsystem, steers?)
# Spin sign follows n_wheel = [0,0,±1] (mirror) in excited_wheel_assembly.dyad;
# spin angle = phi * 1.5 corrects sim tire radius 0.2 -> visual radius 0.4*S.
wheelcfg = [
    (fr_plot, cFR, +1.0, controlled_sys.excited_suspension_fr, true),
    (fl_plot, cFL, -1.0, controlled_sys.excited_suspension_fl, true),
    (br_plot, cBR, +1.0, controlled_sys.excited_suspension_br, false),
    (bl_plot, cBL, -1.0, controlled_sys.excited_suspension_bl, false),
]

function set_pose!(tt)
    T = get_frame(sol, frame_cm, tt)             # 4×4 world transform (NO transpose)
    R = T[1:3, 1:3]
    p = T[1:3, 4]
    tb = OV.Point3f(p + R*o)
    Qb = Qf(R)
    Makie.transform!(body_plot;  translation = tb, rotation = Qb)
    Makie.transform!(glass_plot; translation = tb, rotation = Qb)
    for (plot, ci, dir, susp, steers) in wheelcfg
        phi   = sol(tt, idxs = susp.wheel_rotation.phi)
        spin  = phi * dir * 1.5
        steer = steers ? sol(tt, idxs = susp.steer_rotation.phi) : 0.0
        Rw  = R * Rotations.RotY(steer) * Rotations.RotZ(spin)
        pos = OV.Point3f(p + R*(o .+ ci))
        Makie.transform!(plot; translation = pos, rotation = Qf(Rw))
    end
end

on(t) do tt
    set_pose!(tt)
end
set_pose!(sol.t[1])

# ======================= record (fast no-reset path) ==========================
# OmniverseMakie's offscreen colorbuffer resets RT2 accumulation on any change and
# re-converges `warmup` steps per frame (slow: ~0.4 fps at warmup=32). No-op'ing
# reset! keeps accumulation across frames like the interactive viewport — with
# warmup=4 steps/frame this records ~5-6 fps with no visible motion artifacts.
OV.OV.reset!(r::OV.OV.Renderer; time::Float64 = 0.0) = nothing
logmsg("OV.reset! no-op'd — accumulation carries across frames")

out = joinpath(ROOT, "assets", "controlled_flat_road_car_circle_rt_mesh.mp4")
framerate = 30
timevec = range(sol.t[1], sol.t[end], step = 1/framerate)

for _ in 1:10                                     # pre-roll: warm the t=0 pose
    OV.Makie.colorbuffer(fig)
end
logmsg("pre-roll done; recording $(length(timevec)) frames -> $out")
tstart = time(); nfr = 0
fn = OV.record(fig, out, timevec; framerate) do time
    t[] = time
    global nfr += 1
    nfr % 150 == 0 && logmsg("frame $nfr/$(length(timevec)) ($(round(nfr/(Base.time()-tstart); digits=2)) fps)")
end
logmsg("WROTE $fn size=$(round(filesize(fn)/1e6; digits=1))MB in $(round(time()-tstart; digits=1))s")
