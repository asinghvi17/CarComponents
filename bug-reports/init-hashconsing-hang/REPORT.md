# ODEProblem construction hangs for hours in initialization-problem alias elimination / SymbolicUtils hashconsing

## Summary

Constructing an `ODEProblem` for the closed-loop multibody model
`CarComponents.ControlledFlatRoadCar` (a 4-wheel car with double-wishbone
suspension loop-closures **plus a circular-path steering controller**) hangs on
this machine: a single `ODEProblem(...)` call pegs **one CPU core at ~100% for
7.5+ hours without completing**. The hang is **not** in the ODE right-hand-side
codegen. Captured stacktraces show it is inside the construction of the
**initialization problem** (built eagerly at `ODEProblem` time): building it
triggers `mtkcompile` of the initialization system, whose **alias-elimination**
pass rewrites terms by substitution, and each rebuilt `Add`/`Mul` term is interned
through **SymbolicUtils hashconsing**, whose structural `isequal` on the enormous
multibody expressions explodes (225M+ allocations, only 59 GCs — pure compute in
`isequal`). Bisection (below) pins the trigger to the controller's **algebraic
steering-feedback loop**: the same car plant *without* the feedback (`FullCar`,
`FlatRoadSteeredCar`) builds its `ODEProblem` in ~13-14 s. Passing
`build_initializeprob = false` skips the pathological step entirely (builds in
seconds), confirming the initialization-problem build is the sole locus.

## Environment

- **Julia** 1.12.6, JuliaHub "dyad" distribution (channel `+dyad-3.2.0-next.87`).
- Machine: this JuliaHub workstation. The hang is (at least partly)
  **machine/environment dependent** - the same model reportedly builds in ~19 s on
  a different, fresh machine. It was reproduced **here** (twice, live stacktrace
  captured), where it is catastrophically slow.
- Symbolic / compiler stack - **baked into the sysimage stdlib**
  (`share/julia/stdlib/v1.12/`), i.e. shipped with the distribution, not
  Manifest-resolved:
  - ModelingToolkit 11.29.0, ModelingToolkitBase 1.48.0
  - Symbolics 7.29.0, SymbolicUtils 4.36.0
  - RuntimeGeneratedFunctions 0.5.21
  - SciMLBase 3.24.0
  - DyadCompilerPasses 0.2.6, DyadInterface 7.1.0
- **Manifest-resolved** (loaded from `~/.julia/packages/`):
  - MultibodyComponents 0.2.0 (pinned to tag commit `be76af8`; both `main` and the
    tag reproduce). Provides `World`, `BodyShape`, `Revolute`, `SlippingWheel`,
    `QuarterCarSuspension`, and the `examples.suspension.*` car models used below.
- Sanity control: a trivial 3-equation Lorenz `ODEProblem` builds in **0.5 s** in
  this same environment, so the `ODEProblem` machinery itself is healthy; the
  pathology is specific to this model.

## Symptom

- One `ODEProblem(controlled_sys, defs, (0.0, 30.0); optimize=:basic, saveat=1/30)`
  call runs for **7.5+ hours** without returning.
- **99.8 % of one core**, single-threaded; **stable ~2.58 GB RSS** (not a memory
  blowup / not swapping).
- At interrupt: **`Allocations: 225136742 (Pool: 225128177; Big: 8565); GC: 59`** -
  225 M allocations but only 59 GCs, i.e. almost all time is spent *computing*
  inside `isequal`, not allocating/collecting. Consistent with an O(n^2)-or-worse
  structural comparison. (A second capture, with `chassis_bushings=true`, showed
  the same signature: `Allocations: 236579523; GC: 59`.)
- Contrast: trivial Lorenz `ODEProblem` = **0.5 s**; the *same* model reportedly
  = **~19 s** on a fresh machine.

## Exact reproducer

```julia
using CarComponents, ModelingToolkit, MultibodyComponents, DyadCompilerPasses

@named controlled_car = ControlledFlatRoadCar(
    path_center_x=0.0, path_center_z=-6.0, path_radius=6.0,
    steer_limit=0.6, path_direction=-1.0, heading_gain=1.0, radial_gain=0.15)

reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(;
    inline_linear_sccs=true, analytical_linear_scc_limit=1)

# FAST (~8 s):
controlled_sys = ModelingToolkit.mtkcompile(controlled_car;
    additional_passes=[], reassemble_alg, optimize=[DyadCompilerPasses.LDIV_RULE])

defs = ...  # initial conditions / parameters (see solve_and_save.jl `build_defs`)

# HANGS 7.5+ hours (one core at 100%):
prob = ODEProblem(controlled_sys, defs, (0.0, 30.0); optimize=:basic, saveat=1/30)
```

The full script is `solve_and_save.jl` in this folder. `mtkcompile` of the *main*
system is fast (~8 s, 36 unknowns); the hang is entirely in the subsequent
`ODEProblem` construction.

## Root-cause analysis (from the captured stacktraces)

Two backtrace files accompany this report:

- **`stacktrace_clean.txt`** - read this first; the ~170-frame hashconsing
  `isequal` recursion is collapsed into readable markers.
- **`stacktrace_raw.txt`** - the full 410-line native backtrace.

A **second, independent** live backtrace was captured here when experiment 5
(`chassis_bushings=true`) was SIGTERM'd at the 150 s timeout - it lands in the
identical recursion (`isequal_somescalar` / `isequal_bsimpl` /
`isequal_addmuldict` at `SymbolicUtils/.../hashconsing.jl:37/199/49` under
`find_perfect_aliases!`), corroborating the diagnosis.

The call chain (top-level -> hang), quoting the key frames:

```
ODEProblem                                                   (odeproblem.jl:105)
 -> process_SciMLProblem / __process_SciMLProblem           (problem_utils.jl:1836/1886)
  -> maybe_build_initialization_problem                      (problem_utils.jl:1487/1503)
   -> InitializationProblem                                  (initializationproblem.jl:22)
    -> mtkcompile  (of the INITIALIZATION system)            (systems.jl:84 ...)
     -> _mtkcompile! -> eliminate_perfect_aliases!           (alias_elimination.jl:36)
      -> find_perfect_aliases!                               (alias_elimination.jl:360)
       -> Symbolics Substituter -> IRSubstituter             (equations.jl:37, irstructure.jl:875)
        -> substitute_ir! -> combine_fold -> maketerm        (substitute.jl:222, terminterface.jl:304)
         -> add_worker -> AddWorkerBuffer -> Mul -> AddMul   (addsub.jl:78/128, safe_ctors.jl, inner_ctors.jl:444)
          -> SymbolicUtils hashconsing (WeakCacheSets.getkey!) (WeakCacheSets.jl:201/146)
           -> ht_keyindex2_shorthash! -> isequal -> isequal_addmuldict / isequal_bsimpl
              ...hundreds of frames of structural isequal recursion...   <-- CATASTROPHIC
```

Interpretation:

1. `ODEProblem` construction eagerly builds an **initialization problem**
   (`maybe_build_initialization_problem`).
2. That calls **`mtkcompile` on the initialization system**, which runs
   **alias elimination** (`eliminate_perfect_aliases!` /
   `find_perfect_aliases!`). Alias elimination finds equations of the form
   `a = b` / `a = f(...)` and **substitutes** the RHS everywhere the LHS occurs.
3. Each substituted/rebuilt `Add`/`Mul` term is **interned via SymbolicUtils
   hashconsing** (`WeakCacheSets.getkey!`). Hashconsing must run a structural
   **`isequal`** against cache collisions; on the enormous multibody `Add`/`Mul`
   expressions produced by the analytic loop closures this `isequal` recurses
   hundreds of frames deep and dominates runtime.
4. 225 M+ allocations / 59 GC ⇒ effectively all time is in `isequal` compute. The
   scaling looks super-linear (O(n^2) or worse) in expression size × number of
   alias equations, which is why a modest increase in model coupling (adding the
   controller feedback - see below) moves it from ~13 s to 7.5+ h.

## MWE / bisection results (all reproduced on this machine)

Goal: the smallest self-contained repro. I built up from MultibodyComponents' own
example car models (`using ModelingToolkit, MultibodyComponents`) and bisected
what the controller adds. **Every `ODEProblem` build was run under `timeout 150`**
(a timeout = the hang). Each row is one construction experiment; scripts are
`mwe_attempt_0N_*.jl` in this folder.

| # | Model | Steering joints? | Steer angle source | Controller feedback? | unknowns | `ODEProblem` build | Script |
|---|-------|------------------|--------------------|-----------------------|----------|--------------------|--------|
| 1 | `MultibodyComponents.examples.suspension.FullCar` | no | - | no | 32 | **14.4 s (fast)** | `mwe_attempt_01_fullcar.jl` |
| 2 | `CarComponents.FlatRoadSteeredCar` | yes (front) | **constant parameter** `0.2` | no | 36 | **13.0 s (fast)** | `mwe_attempt_02_steeredcar.jl` |
| 3 | `CarComponents.ControlledFlatRoadCar` | yes (front) | `controller.steer_angle` = f(chassis pose) | **yes** | 36 | **HANG** (mtkcompile 8.2 s, then >80 s in `ODEProblem`, no completion -> 7.5 h) | `mwe_attempt_03_controlledcar.jl` |
| 4 | `ControlledFlatRoadCar`, `build_initializeprob=false` | yes | feedback | yes | 36 | **fast** (skips init build; then errors on missing guesses at +41 s) | `mwe_attempt_04_noinitprob.jl` |
| 5 | `ControlledFlatRoadCar(chassis_bushings=true)` | yes | feedback | yes | 76 | **HANG** (2nd live stacktrace captured; bushings do NOT help) | `mwe_attempt_05_bushings.jl` |

Notes on each:

- **#1 FullCar** is the pure-MultibodyComponents 4-corner plant (4x
  `ExcitedWheelAssembly`, each with a `QuarterCarSuspension` analytic loop
  closure) - i.e. `ControlledFlatRoadCar` *minus the steering + controller*. It
  compiled to 32 unknowns and built its `ODEProblem` in **14.4 s**, emitting the
  init-system "overdetermined (58 eq / 54 unknowns)" warning and a handful of
  Symbolics "Did not converge after maxiters = 100 substitutions" warnings - i.e.
  the alias-substitution machinery is already working hard here, but does **not**
  blow up. **The multibody loop-closure plant alone is not sufficient to trigger
  the hang.**
- **#2 FlatRoadSteeredCar** adds the exact same **steerable** front corners as the
  full model (a `steer_rotation` `Revolute` + a `SteeringPosition` source per
  front wheel), but drives the steering angle from a **constant parameter**. Still
  **fast (13.0 s)**. ⇒ The steerable *structure* is not the trigger either.
- **#3 ControlledFlatRoadCar** differs from #2 by **only** two equations,
  `steer_angle_fr ~ controller.steer_angle` and
  `steer_angle_fl ~ controller.steer_angle`, plus the `CirclePathSteeringController`
  subsystem. The controller is a purely **algebraic** law reading the chassis pose
  (`back_front.body.r_0[1]`, `back_front.body.r_0[3]`, `back_front.body.phi[2]`)
  and producing a steering command through `atan2`, `sqrt`, `atan`, and a `tanh`
  saturation. This closes a **position-level algebraic feedback loop** whose
  right-hand side is a large transcendental expression - exactly the kind of
  `a = f(chassis pose)` perfect-alias equation that alias elimination substitutes
  into the already-huge front-suspension loop terms, tipping hashconsing `isequal`
  into the catastrophic regime. **This is the trigger.**
- **#4** confirms the locus: with `build_initializeprob=false` the constructor
  sails past `maybe_build_initialization_problem` in seconds and instead fails at
  `varmap_to_vars` ("Initial condition underdefined ... missing from the variable
  map" for `...wheel.wheeljoint.delta_0[1]/[3]` and `...rotational_losses.w_rel`).
  So the hang is **entirely** in the initialization-problem build.
- **#5** shows that **decoupling the corners does NOT help**: mounting each corner
  through a compliant 6-DOF `Bushing` (`chassis_bushings=true`) fragments the big
  coupled inline-linear block (76 unknowns), yet the `ODEProblem` build **still
  hangs** in the identical hashconsing recursion. So the blowup is driven by the
  feedback-alias substitution into the (still loop-heavy) initialization system,
  not solely by the single large coupled acceleration block.

**Smallest reproducer found:** `ControlledFlatRoadCar` (the steered car **with**
the controller-feedback equations). The nearest **non**-reproducing neighbour is
`FlatRoadSteeredCar` (identical steerable plant, constant steering) - the two
differ only by the two feedback equations + the algebraic controller. The
pure-MultibodyComponents `FullCar` is the clean negative baseline.

**Not fully minimized (open item):** a *portable* positive reproducer using only
`using ModelingToolkit, MultibodyComponents` would require reconstructing a
steerable corner (a `MultibodyComponents.Revolute` steer joint + position source
inside the `QuarterCarSuspension`/wheel chain) plus a transcendental algebraic
feedback from a body coordinate - essentially re-implementing
`CarComponents.ExcitedWheelAssembly(steering=true)` + `CirclePathSteeringController`
from primitives. MultibodyComponents 0.2.0 does not ship a steerable assembly, so
the minimal *positive* repro here still uses CarComponents' steerable corner. The
decisive, portable finding is the negative baseline + the single-delta A/B:
**multibody loops alone don't hang; adding an algebraic pose->steering feedback loop
does.**

## Suggested workarounds

1. **`ODEProblem(...; build_initializeprob = false)` - VERIFIED to skip the
   hang.** The constructor no longer builds the initialization problem and returns
   in seconds. Caveat (observed): it then needs `u0`/guesses for the variables the
   initialization would otherwise have solved
   (`...wheel.wheeljoint.delta_0[1]/[3]`, `...rotational_losses.w_rel`) or it errors
   with "Initial condition underdefined". So a usable workaround is
   `build_initializeprob = false` **together with** guesses/defaults (or
   `missing_guess_value = MissingGuessValue.Constant(...)`) for those variables.
   *Providing those guesses and confirming a fully-built, solvable problem was not
   completed here - mark UNVERIFIED.*
2. **Relationship to `BrownFullBasicInit()`**: the downstream `solve(...)` in the
   real workflow already passes `initializealg = BrownFullBasicInit()`. Since a
   numerical Brown basic-init is used at solve time anyway, the expensive
   *symbolic* initialization problem built at `ODEProblem` time may be redundant
   for this workflow - making `build_initializeprob = false` (plus guesses) an
   attractive fix. *UNVERIFIED that the solve then succeeds.*
3. **`chassis_bushings = true` - VERIFIED does NOT help** (experiment 5 still
   hangs). Decoupling the corners is not a workaround.
4. **Model-side mitigations still worth trying (UNVERIFIED):**
   - Break the algebraic controller loop by giving the steering command a small
     first-order **filter state** (make `steer_angle` a state driven toward the
     controller output) so the feedback is no longer a position-level *algebraic*
     alias substituted into the loop equations.
   - `wheel_elastic_contact = true` (decouples the wheel-ground contact; different
     mechanism than bushings - untested).
5. **Upstream fixes to consider**: cap/short-circuit structural `isequal` depth in
   SymbolicUtils hashconsing (or trust the hash before a full structural compare);
   and/or make MTK `find_perfect_aliases!` avoid substituting large transcendental
   RHSs into already-large expressions (substitute lazily / bail out above a size
   threshold).

## Where to file

This is fundamentally a **SymbolicUtils hashconsing `isequal` / ModelingToolkit
alias-elimination scaling** problem, surfaced through the initialization-problem
build. The whole symbolic/compiler stack here ships inside the JuliaHub **Dyad**
distribution's sysimage, and the model comes from **MultibodyComponents**
(Manifest-resolved). Recommended:

- File against **JuliaHub DyadIssues** (https://github.com/DyadLang/DyadIssues) -
  the interrupt banner explicitly directs Dyad bug reports there, and the affected
  ModelingToolkit/Symbolics/SymbolicUtils versions are pinned inside the Dyad
  sysimage (so it must be fixed/bumped there). Attach `stacktrace_clean.txt`,
  `stacktrace_raw.txt`, and the `ControlledFlatRoadCar` reproducer.
- Cross-reference upstream **SciML/ModelingToolkit.jl** (alias elimination
  substituting large RHSs) and **JuliaSymbolics/SymbolicUtils.jl** (hashconsing
  `isequal` blow-up), since the root cause lives there.
- Include the A/B bisection (FullCar / FlatRoadSteeredCar fast vs
  ControlledFlatRoadCar hang), the `build_initializeprob=false` data point, and the
  refuted `chassis_bushings=true` mitigation - together they localize the trigger
  precisely.
