# Experiment 04 (WORKAROUND TEST): ControlledFlatRoadCar (the hanging model) but with
# ODEProblem(...; build_initializeprob=false). If this builds fast, it confirms the hang
# is in the initialization-problem construction (alias elim / hashconsing), and gives a workaround.
T0 = time()
using CarComponents, ModelingToolkit, MultibodyComponents, DyadCompilerPasses
logmsg(m) = (println("[+", round(time() - T0; digits = 1), "s] ", m); flush(stdout))
logmsg("packages loaded")

@named controlled_car = ControlledFlatRoadCar(path_center_x=0.0, path_center_z=-6.0, path_radius=6.0, steer_limit=0.6, path_direction=-1.0, heading_gain=1.0, radial_gain=0.15)
logmsg("model built; compiling...")
reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1)
tc = time()
controlled_sys = ModelingToolkit.mtkcompile(controlled_car; additional_passes = [], reassemble_alg, optimize = [DyadCompilerPasses.LDIV_RULE])
logmsg("mtkcompile done in $(round(time() - tc; digits = 1))s  unknowns=$(length(unknowns(controlled_sys)))")

function build_defs(; speed=0.8, init_steer=atan(1/6))
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
    defs
end
defs = build_defs()
logmsg("defs built ($(length(defs)) pairs). Building ODEProblem with build_initializeprob=false ...")
tp = time()
prob = ODEProblem(controlled_sys, defs, (0.0, 30.0); optimize = :basic, saveat = 1/30, build_initializeprob = false)
logmsg("ODEProblem built in $(round(time() - tp; digits = 1))s  u0len=$(length(prob.u0))  [WORKAROUND WORKS]")
logmsg("DONE")
