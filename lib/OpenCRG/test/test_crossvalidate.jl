# lib/OpenCRG/test/test_crossvalidate.jl
using LibOpenCRG

"""
    crossvalidate_grid(path; label=path) -> Int

Load `path` through both this package's `read_crg`/`road_surface_grid` and
the real compiled C oracle (`LibOpenCRG`), then compare the resulting
world-frame grids node-by-node. Returns the number of mismatching nodes
(deliberately NOT a per-node `@test` inside the double loop -- with a real
file's full grid, a systematic bug would otherwise produce hundreds of
near-identical failure messages that bury the actual pattern; the caller
wraps the returned count in a single `@test mismatches == 0`. If that ever
fails, the fix is to drop a `@show i, j, X[i,j], ref_xy.x, Y[i,j], ref_xy.y`
right before the mismatch-counting `if` below to see exactly where and how
it diverges -- not to guess.
"""
function crossvalidate_grid(path::AbstractString)
    data = OpenCRG.read_crg(path)
    u, v, X, Y, Z = OpenCRG.road_surface_grid(data)
    @test size(X) == (length(u), length(v))
    @test size(Y) == size(X)
    @test size(Z) == size(X)

    dsId = crgLoaderReadFile(path)
    @test dsId != 0
    @test crgCheck(dsId) != 0
    cpId = crgContactPointCreate(dsId)
    @test cpId != -1

    mismatches = 0
    for i in eachindex(u), j in eachindex(v)
        ref_xy = crgEvaluv2xy(cpId, u[i], v[j])
        if ref_xy.status != 0
            # NOTE: `≈ ... atol=...` is `@test`-macro sugar, not general Julia
            # syntax -- `!(a ≈ b atol=c)` is a ParseError ("extra tokens after
            # end of expression") when written outside of `@test` (confirmed by
            # direct experiment while implementing this task; the plan's
            # original draft of this loop used the infix form here and would
            # not even parse). `isapprox(a, b; atol=c)` is an ordinary function
            # call and nests inside `!(...)` fine.
            if !isapprox(X[i,j], ref_xy.x; atol=1e-6) || !isapprox(Y[i,j], ref_xy.y; atol=1e-6)
                mismatches += 1
            end
        end
        ref_z = crgEvaluv2z(cpId, u[i], v[j])
        if ref_z.status != 0 && !isnan(ref_z.z) && !isnan(Z[i,j])
            if !isapprox(Z[i,j], ref_z.z; atol=1e-6)
                mismatches += 1
            end
        end
    end

    crgContactPointDelete(cpId)
    crgDataSetRelease(dsId)
    crgMemRelease()
    return mismatches
end

@testset "road_surface_grid, cross-validated against the LibOpenCRG oracle" begin
    # The plan's original fixture list here was just the first two files below
    # -- neither exercises the end-anchored/backward-integration branches of
    # integrate_reference_line/integrate_reference_z (Task 12's review note),
    # and neither exercises lateral_offset_grid's hairpin degeneracy (Task 13's
    # review note). Both gaps are closed by the two synthetic fixtures added
    # below, each independently confirmed to actually load in the real C
    # library (`crgLoaderReadFile` returns a nonzero dsId AND `crgCheck`
    # returns nonzero -- both asserted inside `crossvalidate_grid` above,
    # not just assumed) before any cross-validation result from them is
    # trusted:
    #
    #   - synthetic_end_anchored_2col.crg: like Task 11's synthetic_end_anchored.crg
    #     (end-anchored x/y AND z, non-constant/non-monotonic phi and slope),
    #     but with 3 long-section columns instead of 1 -- the real C library
    #     rejects any file with fewer than 2 (`crgLoader.c`: "no or insufficient
    #     long section data available", checkHeaderConsistency). Cross-validating
    #     this fixture is what actually caught a real, systematic bug in
    #     integrate_reference_z's end-anchored blend (see its docstring in
    #     transform.jl and the updated Task 12 test in test_transform.jl) --
    #     without a >=2-column end-anchored-z fixture, that branch was cross-
    #     validated only against this package's own hand arithmetic, which
    #     reproduced the same (wrong) assumption the implementation made.
    #
    #   - synthetic_hairpin.crg: an exact 180-degree hairpin (equal-length
    #     segments meeting head-on), binary/KDBI so its phi channel encodes
    #     exact `Float64(pi)` with no ASCII-decimal truncation (a first attempt
    #     at this fixture using 7-decimal ASCII/LRFI truncated pi by ~4.6e-8
    #     rad, which turned out to be nowhere near small enough to trigger the
    #     C library's 1.0e-10 epsilon guards -- confirmed by hand before
    #     switching to binary). Cross-validating this fixture confirmed the C
    #     oracle does NOT produce NaN at a hairpin -- it gracefully collapses
    #     to the reference line's own point -- which is what motivated porting
    #     the epsilon guards into lateral_offset_grid (see its docstring) and
    #     updating Task 13's pinning test accordingly.
    #
    #   - belgian_block.crg: the only large, real-world, Float32-binary
    #     fixture in this test suite -- found during the Tasks 16+17 review
    #     to already cross-validate cleanly (0 mismatches) but to have been
    #     missing from this loop the whole time. Added here (unrelated to
    #     the REFPOINT_* fix that's the rest of this change, but cheap to
    #     include in the same pass).
    for fname in ["handmade_curved_minimalist.crg", "handmade_curved_banked_sloped.crg",
                  "synthetic_end_anchored_2col.crg", "synthetic_hairpin.crg", "belgian_block.crg"]
        @testset "$fname" begin
            mismatches = crossvalidate_grid(joinpath(DATA, fname))
            @test mismatches == 0
        end
    end
end
