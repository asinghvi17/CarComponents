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
using GLMakie

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
    prob = @time ODEProblem(controlled_sys, defs, (0.0, tstop); optimize = :basic, saveat = save_dt)
    solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); initializealg = BrownFullBasicInit(), abstol=1e-6, reltol=1e-6, dtmax=0.02)
end

println("Simulating controlled path follower...")
controlled_sol_video = simulate_controlled_car()
println("retcode=", controlled_sol_video.retcode, " samples=", length(controlled_sol_video.t), " t_end=", controlled_sol_video.t[end])

center_x = 0.0; center_z = -6.0; radius = 6.0
rerr(t) = sqrt((controlled_sol_video(t, idxs=controlled_sys.back_front.body.r_0[1])-center_x)^2 + (controlled_sol_video(t, idxs=controlled_sys.back_front.body.r_0[3])-center_z)^2) - radius
errs = [rerr(t) for t in range(5, controlled_sol_video.t[end], length=200)]

out = joinpath((@__DIR__) |> dirname, "assets", "controlled_flat_road_car_circle.mp4")

fn, scene, fig = MultibodyComponents.render(controlled_car, controlled_sol_video;
    filename = out,
    framerate = 30,
    timescale = 1.0,
    x = 8.0,
    y = 12.0,
    z = 8.0,
    lookat = GLMakie.Vec3f(0.0, 0.0, -6.0),
    up = GLMakie.Vec3f(0.0, 1.0, 0.0),
    show_axis = true,
    traces = [controlled_sys.back_front.frame_cm],
    size = (1280, 720),
)