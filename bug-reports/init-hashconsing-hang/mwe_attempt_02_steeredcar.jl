# Experiment 02: FlatRoadSteeredCar = ControlledFlatRoadCar's steerable structure but
# steer_angle driven by a CONSTANT parameter (steer_angle_fr_cmd=0.2) instead of the
# controller feedback. Isolates "steerable structure" from "controller feedback".
# timeout 150 => reproduced hang. Fast build => the controller feedback is the trigger.
T0 = time()
using CarComponents, ModelingToolkit, MultibodyComponents, DyadCompilerPasses
logmsg(m) = (println("[+", round(time() - T0; digits = 1), "s] ", m); flush(stdout))
logmsg("packages loaded")

@named steered_car = FlatRoadSteeredCar()
logmsg("model built; compiling...")
reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1)
tc = time()
sys = ModelingToolkit.mtkcompile(steered_car; additional_passes = [], reassemble_alg, optimize = [DyadCompilerPasses.LDIV_RULE])
logmsg("mtkcompile done in $(round(time() - tc; digits = 1))s  unknowns=$(length(unknowns(sys)))")

function build_defs(; speed=0.8, init_steer=atan(1/6))
    C = [sys.excited_suspension_fr, sys.excited_suspension_fl,
         sys.excited_suspension_br, sys.excited_suspension_bl]
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
        sys.back_front.body.r_0[1] => 0.0,
        sys.back_front.body.r_0[2] => 0.193,
        sys.back_front.body.r_0[3] => 0.0,
        sys.back_front.body.v_0[1] => speed,
        sys.back_front.body.v_0[2] => 0.0,
        sys.back_front.body.v_0[3] => 0.0,
        sys.back_front.body.phi[2] => 0.0,
        sys.excited_suspension_fr.steer_rotation.phi => init_steer,
        sys.excited_suspension_fl.steer_rotation.phi => init_steer,
        sys.excited_suspension_fr.steering_position.w => 0.0,
        sys.excited_suspension_fl.steering_position.w => 0.0,
    ])
    defs
end
defs = build_defs()
logmsg("defs built ($(length(defs)) pairs). Building ODEProblem (optimize=:basic)... [SLOW STEP]")
tp = time()
prob = ODEProblem(sys, defs, (0.0, 30.0); optimize = :basic, saveat = 1/30)
logmsg("ODEProblem built in $(round(time() - tp; digits = 1))s  u0len=$(length(prob.u0))")
logmsg("DONE - no hang")
