# Long-running: build the (pathologically slow here) ODEProblem for ControlledFlatRoadCar,
# solve it, and serialize (model, sys, sol) to disk so the raytracing step can reuse it
# WITHOUT re-solving. Verbose timestamped logging; run detached (nohup) for hours.
using CarComponents, ModelingToolkit, OrdinaryDiffEqDefault, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using ADTypes, SciMLBase, LinearAlgebra, MultibodyComponents, DyadCompilerPasses
using Serialization

const OUT = "/tmp/claude-1000/-home-juliahub-temp-CarComponents/840f304c-875e-4970-bd71-ac38d708413e/scratchpad/car_solution.jls"
logmsg(m) = (println("[", round(time(); digits=1), "] ", m); flush(stdout))

logmsg("imports done")

@named controlled_car = ControlledFlatRoadCar(path_center_x=0.0, path_center_z=-6.0, path_radius=6.0, steer_limit=0.6, path_direction=-1.0, heading_gain=1.0, radial_gain=0.15)
reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1)
t0 = time()
controlled_sys = ModelingToolkit.mtkcompile(controlled_car; additional_passes = [], reassemble_alg, optimize = [DyadCompilerPasses.LDIV_RULE])
logmsg("mtkcompile done in $(round(time()-t0; digits=1))s  unknowns=$(length(unknowns(controlled_sys)))")

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

logmsg("building ODEProblem (optimize=:basic, tstop=30)... [this is the slow step]")
t1 = time()
prob = ODEProblem(controlled_sys, defs, (0.0, 30.0); optimize = :basic, saveat = 1/30)
logmsg("ODEProblem built in $(round(time()-t1; digits=1))s  u0len=$(length(prob.u0))")

# Checkpoint the PROBLEM too (in case solve fails, we keep the expensive build).
try
    serialize(OUT * ".prob", (; controlled_car, controlled_sys, prob))
    logmsg("checkpoint: problem serialized -> $(OUT).prob")
catch e
    logmsg("WARN: prob serialize failed: $e")
end

logmsg("solving (Rodas5P, 0-30s)...")
t2 = time()
sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); initializealg = BrownFullBasicInit(), abstol=1e-6, reltol=1e-6, dtmax=0.02)
logmsg("solve done in $(round(time()-t2; digits=1))s  retcode=$(sol.retcode)  samples=$(length(sol.t))  t_end=$(sol.t[end])")

logmsg("serializing solution...")
try
    serialize(OUT, (; controlled_car, controlled_sys, sol))
    logmsg("SERIALIZED OK -> $OUT  ($(round(filesize(OUT)/1e6; digits=1)) MB)")
catch e
    logmsg("ERROR: solution serialize failed: $e")
    # fallback: raw arrays
    try
        serialize(OUT * ".raw", (; t = sol.t, u = sol.u, retcode = sol.retcode))
        logmsg("fallback raw arrays serialized -> $(OUT).raw")
    catch e2
        logmsg("ERROR: raw serialize also failed: $e2")
    end
end
logmsg("ALL DONE")
