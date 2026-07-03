# MWE attempt 01: FullCar (ControlledFlatRoadCar minus the steering controller).
# Same compile path: multibody() = mtkcompile(analytical_linear_scc_limit=1, LDIV_RULE),
# then ODEProblem(...; optimize=:basic). Only needs ModelingToolkit + MultibodyComponents.
# Wrap in `timeout 150`; timeout during ODEProblem == reproduced the init/hashconsing hang.
T0 = time()
using ModelingToolkit, MultibodyComponents
logmsg(m) = (println("[+", round(time() - T0; digits = 1), "s] ", m); flush(stdout))
logmsg("packages loaded")

const susp = MultibodyComponents.examples.suspension

logmsg("building FullCar model...")
@named fc = susp.FullCar()
logmsg("model built; compiling (multibody)...")
tc = time()
fc_sys = multibody(fc)
logmsg("mtkcompile done in $(round(time() - tc; digits = 1))s  unknowns=$(length(unknowns(fc_sys)))")

C = [fc_sys.excited_suspension_fr, fc_sys.excited_suspension_fl,
     fc_sys.excited_suspension_br, fc_sys.excited_suspension_bl]
defs = Pair[]
for c in C
    push!(defs, c.wheel.wheeljoint.v_small => 1e-3)
    push!(defs, c.suspension.ks => 5 * 44000)
    push!(defs, c.suspension.cs => 5 * 4000)
    push!(defs, c.suspension.r2.phi => 5.932380614359173)
    push!(defs, c.wheel_rotation.phi => 0.0)
    push!(defs, c.wheel_rotation.w => 0.0)
end
append!(defs, [
    fc_sys.back_front.body.r_0[1] => 0.0
    fc_sys.back_front.body.r_0[2] => 0.193
    fc_sys.back_front.body.r_0[3] => 0.0
    fc_sys.back_front.body.v_0[1] => 1.0
])
logmsg("defs built ($(length(defs)) pairs). Building ODEProblem (optimize=:basic)... [SLOW STEP]")
tp = time()
prob = ODEProblem(fc_sys, defs, (0, 1.0); optimize = :basic)
logmsg("ODEProblem built in $(round(time() - tp; digits = 1))s  u0len=$(length(prob.u0))")
logmsg("DONE - no hang")
