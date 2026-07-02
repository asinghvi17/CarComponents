# CarComponents

CarComponents is a Dyad / ModelingToolkit workspace for multibody vehicle models.
It currently focuses on a four-wheel car with suspension, slipping tire contact,
flat-road excitation, steering inputs, and a separate continuous-time path-following
steering controller.

The generated Julia files in `generated/` are compiler output. Do not edit them
directly; edit the `.dyad` source files in `dyad/` and recompile.

## Model structure

The Dyad sources are split by responsibility:

| File | Main definitions |
| --- | --- |
| `dyad/slipping_wheel.dyad` | `SlippingWheel`, a local wrapper around `MultibodyComponents.SlipWheelJoint` |
| `dyad/steering_position.dyad` | `SteeringPosition`, a filtered rotational steering-position source |
| `dyad/excited_wheel_assembly.dyad` | `ExcitedWheelAssembly`, suspension + wheel + road-height input |
| `dyad/full_car.dyad` | `FullCar`, the four-corner vehicle plant with external road-height and steering ports |
| `dyad/flat_road_steered_car.dyad` | `FlatRoadSteeredCar`, flat road with fixed front steering commands |
| `dyad/path_controller.dyad` | `CirclePathSteeringController` and `ControlledFlatRoadCar` |
| `dyad/hello.dyad` | Basic starter example |

## Main components

### `FullCar`

`FullCar` is the plant model. It keeps the multibody `World` at the top level and
exposes the intended causal ports:

- road-height inputs for all four wheels:
  - `road_height_fr`
  - `road_height_fl`
  - `road_height_br`
  - `road_height_bl`
- front steering inputs:
  - `steer_angle_fr`
  - `steer_angle_fl`
- wheel contact position outputs:
  - `wheel_position_fr`
  - `wheel_position_fl`
  - `wheel_position_br`
  - `wheel_position_bl`

The vertical world coordinate is `Y`; the road/path plane is `X-Z`.

### `FlatRoadSteeredCar`

`FlatRoadSteeredCar` extends `FullCar`, fixes all road heights to zero, and drives
the two front steering inputs from parameters. This is useful for open-loop
steering demonstrations.

### `CirclePathSteeringController`

`CirclePathSteeringController` is a causal controller separated from the plant. It
reads vehicle position and heading in the world `X-Z` plane and produces a front
steering command.

Inputs:

- `x`
- `z`
- `heading`

Output:

- `steer_angle`

Path parameters include:

- `center_x`
- `center_z`
- `radius`
- `direction`
- `heading_gain`
- `radial_gain`
- `steer_limit`

The controller uses tangent-heading feedforward plus heading and radial-error
feedback with smooth steering saturation.

### `ControlledFlatRoadCar`

`ControlledFlatRoadCar` extends `FullCar` and instantiates a separate
`CirclePathSteeringController`:

```text
FullCar plant + flat road  <->  CirclePathSteeringController
```

The wrapper connects chassis pose to the controller and the controller output to
both front steering inputs:

```dyad
controller.x = back_front.body.r_0[1]
controller.z = back_front.body.r_0[3]
controller.heading = -back_front.body.phi[2]

steer_angle_fr = controller.steer_angle
steer_angle_fl = controller.steer_angle
```

## Typical Julia workflow

Load the package and compile a multibody model with multibody-oriented structural
settings:

```julia
using CarComponents
using ModelingToolkit
using MultibodyComponents
using DyadCompilerPasses

@named car = ControlledFlatRoadCar()

reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(;
    inline_linear_sccs = true,
    analytical_linear_scc_limit = 1,
)

sys = ModelingToolkit.mtkcompile(
    car;
    additional_passes = [],
    reassemble_alg,
    optimize = [DyadCompilerPasses.LDIV_RULE],
)
```

Then build and solve an `ODEProblem` with representative initial conditions for
chassis speed, suspension states, and wheel spin.

## Rendering

Multibody rendering uses `MultibodyComponents.render` and a Makie backend:

```julia
using GLMakie

MultibodyComponents.render(
    car, sol;
    filename = "assets/controlled_flat_road_car_circle.mp4",
    framerate = 30,
    show_axis = true,
    traces = [sys.back_front.frame_cm],
)
```

Recent rendered examples are stored under `assets/`, including videos of:

- open-loop flat-road steering
- controlled circular-path tracking
- controlled path tracking with desired/actual path overlays

## Road-data interpolation direction

The workspace includes `DataInterpolationsND` as a dependency for planned ND road
and map interpolation. The intended integration pattern is:

1. Load or construct grid data in Julia.
2. Build a `DataInterpolationsND.NDInterpolation` object in a Julia helper.
3. Pass that interpolation object into a Dyad component as a native parameter.
4. Evaluate it inside relations, e.g. `height = table(x, z)`.

A future road subsystem can then keep road data separate from the car plant:

```text
FullCar wheel X/Z positions  ->  interpolated road sampler  ->  road-height inputs
```

`lib/OpenCRG` is a work-in-progress reader for OpenCRG road files. The car models
currently use flat road inputs and do not depend on OpenCRG data at runtime.

## Development notes

- Edit `.dyad` files in `dyad/`.
- Recompile after Dyad edits.
- Treat `generated/` as read-only compiler output.
- Keep plant, controller, and road-data functionality in separate components so
they can be tested independently.
