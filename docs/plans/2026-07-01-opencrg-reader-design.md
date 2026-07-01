# OpenCRG Reader — Design

## Context

CarComponents currently drives per-wheel road height (`ExcitedWheelAssembly.road_height`) from a synthetic sine-bump generator (`RoadData` / `dyad/full_car_with_road_ports.dyad`). This design adds a Julia reader for [ASAM OpenCRG](https://github.com/asam-ev/OpenCRG) `.crg` road-surface-measurement files, so real/recorded road surfaces can eventually drive that same port. `Project.toml` already carries a placeholder `lib/OpenCRG` source dependency (currently a "Hello World" stub) and an existing `DataInterpolationsND` dependency that this design builds toward using.

## Two-package hybrid strategy

- **`lib/OpenCRG`** (this design's subject): a pure-Julia, stdlib-only parser + forward-transform. This is the package CarComponents will actually depend on.
- **`lib/LibOpenCRG`** (already implemented): a `ccall` wrapper around the upstream C reference library (github.com/asam-ev/OpenCRG, Apache-2.0), vendored source (`csrc/`) compiled locally via a flat `cc -shared` invocation into a gitignored `build/` directory (rebuild via `julia lib/LibOpenCRG/build.jl`). Its only role is as a correctness oracle in `lib/OpenCRG`'s test suite — it is not meant for production use in vehicle models. It binds 23 functions from `crgBaseLib.h`, including the loader, contact-point evaluation (`crgEvaluv2xy`/`crgEvalxy2uv`/`crgEvaluv2z`/`crgEvalxy2z`/`crgEvaluv2pk`/`crgEvalxy2pk`), and range/increment queries. Its own smoke test (65 assertions) already verifies `crgEvaluv2xy`/`crgEvalxy2uv` round-trip to 1e-6 on the vendored example files.

## File format background

An OpenCRG file is a plain-text header of `$SECTION` blocks (`$CT`, `$ROAD_CRG`, `$ROAD_CRG_OPTS`, `$ROAD_CRG_MODS`, `$ROAD_CRG_MPRO`, `$ROAD_CRG_FILE`, `$KD_DEFINITION`), terminated by a `$$$$...` line, followed by a data grid. The grid is a road cross-section (`v`, lateral offset) repeated at every longitudinal step (`u`, from `REFERENCE_LINE_INCREMENT`), stored in one of four encodings (`LRFI`/`LDFI`/`KRBI`/`KDBI` — ASCII/binary × float32/float64), each a fixed-width/fixed-record-size layout. Per-row channels (`phi` heading, `banking`, `slope`) ride alongside the elevation grid. Full details in the upstream `doc/specification/` (AsciiDoc/Antora) tree.

## What `lib/OpenCRG` parses

The complete header (all sections listed above) and all four data-payload encodings.

## What `lib/OpenCRG` computes

A single batched, vectorized forward transform, `road_surface_grid`:
1. Integrate heading (`phi`) via `cumsum` along `u` to get the reference line's world-frame path `x_ref(u)`, `y_ref(u)` (accounting for the spec's row-shift: row *i*'s `phi` is the heading of the segment arriving at point *i*, row 0 is a placeholder).
2. Broadcast each `v` offset perpendicular to local heading to get full `X[iu,iv]`, `Y[iu,iv]` grids.
3. Combine the `z(u,v)` elevation grid with banking/slope contributions into `Z[iu,iv]`.
4. Apply `$ROAD_CRG_MODS` (scale/offset/rotate/translate) as part of this same transform.

This is forward-only — no per-point Newton-iteration inversion (`(x,y) → (u,v)`, the numerically fiddly part of OpenCRG) is implemented, since the whole grid is transformed at once rather than queried point-by-point.

## Parsed but inert (explicit non-goals for v1)

- **`$ROAD_CRG_MPRO`** (ellipsoid/datum/Helmert-transform/projection metadata — a geospatial CRS, not a matrix): stored verbatim, no projection math implemented. This only matters for placing the road on a real-world map, which is disconnected from computing road height under a wheel.
- **`$ROAD_CRG_OPTS`** (border modes, warning levels, check tolerances): stored; border/extrapolation behavior is left to the consumer (e.g. `DataInterpolationsND.ExtrapolationType`), not implemented here.
- No interpolation/evaluation API at all — `lib/OpenCRG` returns grids, not an evaluator.

## Package boundary

`lib/OpenCRG` has no dependency on `DataInterpolationsND`, `DyadData`, or `ModelingToolkit` — it returns plain Julia arrays. Wiring `(u, v, X, Y, Z)` into an `NDInterpolation` (two `LinearInterpolationDimension`s, vector-valued `[X,Y,Z]` output) or a Dyad component is CarComponents' job, kept separate so `OpenCRG.jl` stays standalone, independently testable, and potentially publishable on its own.

## Public API (sketch)

```julia
read_crg(path) -> CRGData
road_surface_grid(data::CRGData) -> (u, v, X, Y, Z)
```

`CRGData` holds the parsed reference-line parameters, raw channels (`phi`, `banking`, `slope`, `z`), and a catch-all store for `$ROAD_CRG_OPTS`/`$ROAD_CRG_MPRO` fields.

## Testing strategy

- Vendor the small Apache-2.0-licensed example files already sitting in `lib/LibOpenCRG/test/data/` (`handmade_curved_minimalist.crg`, `handmade_curved_banked_sloped.crg`) for round-trip parsing tests.
- Cross-validate `road_surface_grid`'s `(X,Y)` output against `LibOpenCRG.crgEvaluv2xy` at sample `(u,v)` points.
- Round-trip test each of the 4 payload encodings against the same logical grid (if example files covering all 4 aren't available upstream, synthesize small ones from a parsed ASCII file).
