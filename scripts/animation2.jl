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

@named controlled_car = ControlledFlatRoadCar(path_center_x=0.0, path_center_z=-6.0, path_radius=6.0, steer_limit=0.6, path_direction=-1.0, heading_gain=1.0, radial_gain=0.15)
reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1)
println("Compiling controlled car...")
controlled_sys = ModelingToolkit.mtkcompile(controlled_car; additional_passes = [], reassemble_alg, optimize = [DyadCompilerPasses.LDIV_RULE])
println("compiled unknowns=", length(unknowns(controlled_sys)), " equations=", length(equations(controlled_sys)))

function simulate_controlled_car(; speed=0.8, tstop=30.0, save_dt=1/30, init_steer=atan(1/6))
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
        # Seed the algebraic unknowns that the symbolic initialization problem would
        # otherwise solve for. Required because we build with build_initializeprob=false
        # below (see NOTE); BrownFullBasicInit refines these to consistency at solve time.
        push!(defs, c.wheel.wheeljoint.delta_0[1] => 0.0)
        push!(defs, c.wheel.wheeljoint.delta_0[3] => 0.0)
        push!(defs, c.rotational_losses.w_rel => 0.0)
        # There is no drive torque in the model — the car coasts and tire-slip scrub
        # (not this damper) parks it by ~t=35-40s. Reduced losses kept for the tail.
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
    # NOTE: build_initializeprob=false is REQUIRED on this machine. The default
    # (symbolic) initialization-problem build triggers a SymbolicUtils hashconsing
    # blowup (mtkcompile of the init system -> alias elimination -> isequal_addmuldict
    # recursion) that hangs ODEProblem construction for 7.5+ hours here (though it is
    # ~instant on some machines). Skipping it drops construction to ~4s; the 12 algebraic
    # unknowns seeded above are then made consistent numerically by BrownFullBasicInit.
    # See bug-reports/init-hashconsing-hang/ for the full diagnosis.
    prob = ODEProblem(controlled_sys, defs, (0.0, tstop); optimize = :basic, saveat = save_dt, build_initializeprob = false)
    solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); initializealg = BrownFullBasicInit(), abstol=1e-6, reltol=1e-6, dtmax=0.02)
end

println("Simulating controlled path follower...")
# 45s: the car coasts to rest ON the circle by ~t=40 (no drive torque in the model);
# a longer tstop just films a parked car.
controlled_sol_video = simulate_controlled_car(tstop = 45.0)
println("retcode=", controlled_sol_video.retcode, " samples=", length(controlled_sol_video.t), " t_end=", controlled_sol_video.t[end])

center_x = 0.0; center_z = -6.0; radius = 6.0
rerr(t) = sqrt((controlled_sol_video(t, idxs=controlled_sys.back_front.body.r_0[1])-center_x)^2 + (controlled_sol_video(t, idxs=controlled_sys.back_front.body.r_0[3])-center_z)^2) - radius
errs = [rerr(t) for t in range(5, controlled_sol_video.t[end], length=200)]

out = joinpath((@__DIR__) |> dirname, "assets", "controlled_flat_road_car_circle_rt.mp4")

# --- Raytraced rendering via OmniverseMakie ---------------------------------
# warmup = RTX accumulation steps per frame (quality/speed knob; `samples` is a no-op
# for offscreen rendering as of OmniverseMakie 0.1.0 — see ARCHITECTURE.md).
OmniverseMakie.activate!(warmup = 32)

# A flat ground quad so the car has something to sit on — the sim has no visible
# road geometry of its own (road_height is just a scalar input, always 0 here).
ground_half = radius + 8.0
ground = OmniverseMakie.Rect3f(
    OmniverseMakie.Vec3f(center_x - ground_half, -0.05, center_z - ground_half),
    OmniverseMakie.Vec3f(2ground_half, 0.05, 2ground_half),
)

# 3-point-ish lighting: dim cool ambient dome (no HDR env maps in OmniverseMakie yet)
# + warm key + cool fill + subtle rim. PointLight is (color, position), color-first.
# Intensities deliberately low — hotter values blow out the white floor.
env_bg = fill(OmniverseMakie.RGBf(0.30, 0.33, 0.38), 1, 1)
lights = [
    OmniverseMakie.EnvironmentLight(0.22, env_bg'),
    OmniverseMakie.PointLight(OmniverseMakie.RGBf(240, 205, 165), OmniverseMakie.Vec3f(center_x + 5, 16, center_z + 4)),
    OmniverseMakie.PointLight(OmniverseMakie.RGBf(55, 70, 100), OmniverseMakie.Vec3f(center_x - 9, 8, center_z - 12)),
    OmniverseMakie.PointLight(OmniverseMakie.RGBf(95, 85, 105), OmniverseMakie.Vec3f(center_x - 2, 10, center_z + 12)),
]

# MultibodyComponents.render's animation entry (render(model, sol; filename, ...)) opens
# its own record loop internally with no hook to add extra geometry first. Use the
# still-frame entry instead (render(model, sol, time::Real; ...)) which returns the
# live (fig, t, scene) before recording, add the ground quad, then drive `t` ourselves
# in the same way the animation entry does internally (ext/Render.jl:763).
fig, t, scene = MultibodyComponents.render(controlled_car, controlled_sol_video, controlled_sol_video.t[1];
    slider = false,
    x = 4.2,   # zoomed in vs animation1's (8, 12, 8)
    y = 6.2,
    z = 4.2,
    lookat = OmniverseMakie.Vec3f(0.0, 0.0, -6.0),
    up = OmniverseMakie.Vec3f(0.0, 1.0, 0.0),
    traces = [controlled_sys.back_front.frame_cm],
    size = (1280, 720),
    lights = lights,
)
ground_plot = OmniverseMakie.mesh!(scene, ground; color = :gainsboro)

# Reference (target) circle the controller tracks, as a thin floor-level ring —
# the cyan trace is the REAL path; this green ring is the intended one.
let θs = range(0, 2π, length = 260)
    ref_pts = [OmniverseMakie.Point3f(center_x + radius * cos(a), 0.03, center_z + radius * sin(a)) for a in θs]
    OmniverseMakie.lines!(scene, ref_pts; color = OmniverseMakie.RGBf(0.12, 0.8, 0.2), linewidth = 1.2)
end

# --- Per-part PBR materials (post-hoc on the Render.jl-created Mesh plots) ---------
# Render.jl draws matte displayColor primitives. OmniverseMakie tracks plot.material[]
# live, so classify every Mesh (except the ground) by its color and assign OmniPBR:
# pure-red tire discs -> dark rubber, light-grey rims -> brushed metal, everything
# else -> red body paint. Paint is DIELECTRIC (metallic 0): metallic paint reflects
# the bright floor and washes out to pink.
let
    meshes = Any[]
    stack = Any[scene.scene]
    while !isempty(stack)
        s = pop!(stack)
        for p in s.plots
            (p === ground_plot) && continue
            p isa OmniverseMakie.Mesh && push!(meshes, p)
        end
        append!(stack, s.children)
    end
    for p in meshes
        c = try
            OmniverseMakie.to_color(p.color[])
        catch
            OmniverseMakie.RGBAf(0.7, 0.7, 0.7, 1)
        end
        r, g, b = Float32(OmniverseMakie.red(c)), Float32(OmniverseMakie.green(c)), Float32(OmniverseMakie.blue(c))
        if r > 0.55 && g < 0.45 && b < 0.45
            p.material[] = (; base_color = OmniverseMakie.RGBf(0.03, 0.03, 0.035), metallic = 0.0f0, roughness = 0.85f0)
        elseif min(r, g, b) > 0.5
            p.material[] = (; metallic = 0.9f0, roughness = 0.25f0)
        else
            p.material[] = (; base_color = OmniverseMakie.RGBf(0.45, 0.02, 0.03), metallic = 0.0f0, roughness = 0.28f0)
        end
    end
end

framerate = 30
timescale = 1.0
timevec = range(controlled_sol_video.t[1], controlled_sol_video.t[end] * timescale, step = 1 / framerate)
println("Recording ", length(timevec), " frames to ", out)
fn = OmniverseMakie.record(fig, out, timevec; framerate) do time
    t[] = time / timescale
end
println("Wrote ", fn)
