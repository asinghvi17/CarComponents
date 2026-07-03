# Experiment 05 (WORKAROUND/localization): ControlledFlatRoadCar with chassis_bushings=true.
# Per the component docs, this mounts each corner through a compliant 6-DOF Bushing, decoupling
# the four corners so the single large coupled inline-linear block fragments into per-corner blocks.
# If ODEProblem then builds fast, it confirms the hashconsing blowup is driven by the big coupled
# block being fed the transcendental steering-feedback substitution.
T0 = time()
using CarComponents, ModelingToolkit, MultibodyComponents, DyadCompilerPasses
logmsg(m) = (println("[+", round(time() - T0; digits = 1), "s] ", m); flush(stdout))
logmsg("packages loaded")

@named controlled_car = ControlledFlatRoadCar(path_center_x=0.0, path_center_z=-6.0, path_radius=6.0,
    steer_limit=0.6, path_direction=-1.0, heading_gain=1.0, radial_gain=0.15,
    chassis_bushings=true)
logmsg("model built (chassis_bushings=true); compiling...")
reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1)
tc = time()
controlled_sys = ModelingToolkit.mtkcompile(controlled_car; additional_passes = [], reassemble_alg, optimize = [DyadCompilerPasses.LDIV_RULE])
logmsg("mtkcompile done in $(round(time() - tc; digits = 1))s  unknowns=$(length(unknowns(controlled_sys)))")

# Minimal defs: chassis pose + per-corner core; rely on missing_guess_value for the rest.
defs = Pair[
    controlled_sys.back_front.body.r_0[1] => 0.0,
    controlled_sys.back_front.body.r_0[2] => 0.193,
    controlled_sys.back_front.body.r_0[3] => 0.0,
    controlled_sys.back_front.body.v_0[1] => 0.8,
    controlled_sys.back_front.body.phi[2] => 0.0,
]
logmsg("Building ODEProblem (optimize=:basic, missing_guess_value=Constant(0.0))...")
tp = time()
prob = ODEProblem(controlled_sys, defs, (0.0, 30.0); optimize = :basic, saveat = 1/30,
    missing_guess_value = MissingGuessValue.Constant(0.0))
logmsg("ODEProblem built in $(round(time() - tp; digits = 1))s  u0len=$(length(prob.u0))  [NO HANG with bushings]")
logmsg("DONE")
