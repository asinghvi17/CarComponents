# lib/OpenCRG/test/test_transform.jl
@testset "integrate_reference_line" begin
    @testset "Case A: no end anchoring, simple forward Euler" begin
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_START_X=0.0", "REFERENCE_LINE_START_Y=0.0"])
        phi = [0.0, 0.0, 0.0, 0.0]   # straight line along +x; phi[1] unused
        x, y = OpenCRG.integrate_reference_line(r, phi)
        @test x ≈ [0.0, 1.0, 2.0, 3.0]
        @test y ≈ [0.0, 0.0, 0.0, 0.0]
    end

    @testset "Case B: end-anchored, error redistributed linearly" begin
        data = OpenCRG.read_crg(joinpath(DATA, "synthetic_end_anchored.crg"))
        x, y = OpenCRG.integrate_reference_line(data.refline, data.phi)
        @test x[1] ≈ 0.0
        @test x[end] ≈ 4.4          # true end position is hit exactly
        # 0.4m of error over 4 segments, redistributed with fraction i/(n-1).
        # Backward pass (phi≡0, du=1): xb = [end_x - k for k in 4:-1:0] = [0.4, 1.4, 2.4, 3.4, 4.4].
        # Forward/blend recursion (NOTE: chains off the already-blended x[i], not
        # off a separate unblended forward trajectory, so the correction compounds
        # instead of landing on the "obviously linear" 0, 1.1, 2.2, 3.3, 4.4):
        #   i=1, frac=1/4: x[2] = 0.75*(x[1]+1) + 0.25*xb[2] = 0.75*1 + 0.25*1.4 = 1.1
        #   i=2, frac=2/4: x[3] = 0.5*(x[2]+1) + 0.5*xb[3]  = 0.5*2.1 + 0.5*2.4 = 1.05 + 1.2 = 2.25
        # The plan's original hand-computed literal here (2.1333333333333333) was
        # independently checked -- by hand, in an isolated Julia scratch run, and
        # via this exact end-to-end read_crg/integrate_reference_line pipeline --
        # and all three agree the correct value is 2.25, not 2.1333333333333333.
        @test x ≈ [0.0, 1.1, 2.25, 3.3625, 4.4] atol=1e-12   # full vector, not just node 3 -- x[2]/x[4] compound into x[3], but assert them all directly too
        @test all(y .≈ 0.0)
    end

    @testset "arrival-heading convention is the same in both the forward and backward passes" begin
        # phi≡0 everywhere structurally cannot distinguish phi[i] from phi[i+1]
        # (cos(0)=1, sin(0)=0 regardless of index), so a backward-pass index
        # bug would slip through the tests above undetected. Use a genuinely
        # curving phi instead: integrate it forward-only (Case A), then re-run
        # Case B anchored to Case A's own true endpoint -- if the backward
        # pass reads the wrong index, this will NOT reproduce Case A exactly.
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0"])
        phi = [0.0, 0.1, 0.3, 0.2, -0.1, 0.4]   # phi[1] unused; deliberately non-constant and non-monotonic
        xa, ya = OpenCRG.integrate_reference_line(r, phi)

        r_anchored = OpenCRG.ReferenceLineParams(
            r.start_u, r.end_u, r.increment, r.start_x, r.start_y, r.start_phi,
            xa[end], ya[end], r.end_phi, r.start_z, r.end_z,
            r.v_right, r.v_left, r.v_increment, r.start_slope, r.end_slope, r.start_banking, r.end_banking,
        )
        xb, yb = OpenCRG.integrate_reference_line(r_anchored, phi)
        @test xa ≈ xb atol=1e-12
        @test ya ≈ yb atol=1e-12
    end

    @testset "exactly one of end_x/end_y set is an error" begin
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_END_X=4.4"])
        @test r.end_y === nothing   # sanity: END_Y really wasn't declared
        @test_throws Exception OpenCRG.integrate_reference_line(r, [0.0, 0.0, 0.0])
    end
end

@testset "integrate_reference_z" begin
    @testset "no slope channel, zero start_z: early-out to all zeros" begin
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0"])
        @test OpenCRG.integrate_reference_z(r, nothing, 4) == zeros(4)
    end

    @testset "no slope channel, NONZERO start_z: early-out to constant start_z, not zero" begin
        # Caught by code review, cross-checked against the real compiled C
        # reference library (crgEvalz.c's fallback: crgData->channelRefZ.info.first,
        # i.e. REFERENCE_LINE_START_Z, not zero, when the ref-z channel is invalid).
        # A flat road at a nonzero elevation must report that elevation everywhere,
        # not silently flatten to 0.0.
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_START_Z=5.0"])
        @test OpenCRG.integrate_reference_z(r, nothing, 4) == fill(5.0, 4)
    end

    @testset "constant slope (no channel, nonzero REFERENCE_LINE_START_S)" begin
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_START_S=0.1"])
        z_ref = OpenCRG.integrate_reference_z(r, nothing, 4)
        @test z_ref ≈ [0.0, 0.1, 0.2, 0.3]
    end

    @testset "per-row slope channel, no end anchoring" begin
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_START_Z=1.0"])
        slope = [0.0, 0.1, 0.2, 0.3]   # slope[1] unused, matching the phi convention
        z_ref = OpenCRG.integrate_reference_z(r, slope, 4)
        @test z_ref ≈ [1.0, 1.1, 1.3, 1.6]
    end

    @testset "Case B: end-anchored, non-constant slope so the blend is actually checkable" begin
        # The plan as originally written has NO test at all for the
        # refline.end_z !== nothing (backward-integrate-then-blend) branch. Worse,
        # the most "natural" fixture to bolt on by analogy with Task 11 would reuse
        # a constant slope (or slope===nothing with nonzero start_slope) -- but that
        # is exactly Task 11's phi≡0 blind spot in a new costume: slope_at(i) and
        # slope_at(i-1) (or any other off-by-one) return the SAME constant
        # regardless of index, so a backward-pass index bug is structurally
        # undetectable. Confirmed by injecting a `slope_at(i-1)` typo (should be
        # `slope_at(i)`) into a scratch copy of `integrate_reference_z`: with
        # constant slope, correct and buggy outputs were bit-for-bit identical;
        # with the non-constant slope used below, they diverge sharply (e.g.
        # z_ref[2] 0.2 correct vs 0.25 buggy). Hence: non-constant slope here.
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_START_Z=0.0", "REFERENCE_LINE_END_Z=0.9"])
        slope = [NaN, 0.1, 0.2, -0.1, 0.3]   # slope[1]=NaN: never read: a stray read would poison the whole result
        z_ref = OpenCRG.integrate_reference_z(r, slope, 5)
        @test z_ref[1] ≈ 0.0
        @test z_ref[end] ≈ 0.9   # true end elevation is hit exactly, same as Task 11's Case B for (x,y)
        # Backward pass (du=1, unaffected by the Task 17 fix below): zb[5]=0.9,
        # zb[4]=0.9-0.3=0.6, zb[3]=0.6-(-0.1)=0.7, zb[2]=0.7-0.2=0.5,
        # zb[1]=0.5-0.1=0.4 -> zb=[0.4, 0.5, 0.7, 0.6, 0.9].
        #
        # UPDATED during Task 17's cross-validation: the literals below were
        # originally [0.0, 0.2, 0.55, 0.5625, 0.9], derived (both by hand and by
        # this package's own from-scratch implementation) by exact analogy with
        # Task 11's (x,y) blend formula, `fraction = i/(nu-1)`. Widening Task 17's
        # cross-validation fixture set to include a >=2-long-section-column
        # end-anchored file (the original 1-column synthetic_end_anchored.crg
        # can't even load in the real C library) surfaced a real, systematic
        # ~0.1-0.15m mismatch against the compiled oracle at every interior node
        # -- traced to `calcRefLineZ` (crgLoader.c) NOT being symmetric with
        # `calcRefLine`: it uses `fraction = i/(size-1)` (the raw 0-based loop
        # index, not `i+1`) and stops one iteration earlier, so the true
        # `fraction` at a given node is one node "behind" what the (x,y) case
        # uses. See `integrate_reference_z`'s docstring in transform.jl for the
        # full derivation. Corrected recursion (fraction now uses the PRIOR
        # node's position in the blend, i.e. `(i-1)/(nu-1)` for the loop's `i`):
        #   z_ref[5] = zb[5] = 0.9                              (anchor, set directly, never re-blended)
        #   i=1, frac=0/4=0:    z_ref[2] = 1.00*(0.0+0.1)  + 0.00*0.5 = 0.1
        #   i=2, frac=1/4=0.25: z_ref[3] = 0.75*(0.1+0.2)  + 0.25*0.7 = 0.225 + 0.175 = 0.4
        #   i=3, frac=2/4=0.5:  z_ref[4] = 0.50*(0.4-0.1)  + 0.50*0.6 = 0.15  + 0.3   = 0.45
        # Independently verified three ways: by hand (above), in an isolated
        # Julia scratch session, and by cross-validating against the compiled C
        # oracle's actual `crgEvaluv2z` output on `synthetic_end_anchored_2col.crg`
        # (see test_crossvalidate.jl), which agrees bit-for-bit (within 1e-6).
        @test z_ref ≈ [0.0, 0.1, 0.4, 0.45, 0.9] atol=1e-12
    end

    @testset "arrival-slope convention is the same in both the forward and backward passes" begin
        # Same round-trip idea as Task 11's analogous phi test: integrate a
        # genuinely non-constant, non-monotonic slope forward-only (no end_z) to
        # get a true endpoint, then re-run anchored to that exact endpoint -- if
        # the backward pass reads the wrong slope index, this will NOT reproduce
        # the forward-only result.
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_START_Z=2.0"])
        slope = [NaN, 0.05, -0.2, 0.15, 0.4, -0.05]   # slope[1] unused; non-constant and non-monotonic
        za = OpenCRG.integrate_reference_z(r, slope, 6)

        r_anchored = OpenCRG.ReferenceLineParams(
            r.start_u, r.end_u, r.increment, r.start_x, r.start_y, r.start_phi,
            r.end_x, r.end_y, r.end_phi, r.start_z, za[end],
            r.v_right, r.v_left, r.v_increment, r.start_slope, r.end_slope, r.start_banking, r.end_banking,
        )
        zb = OpenCRG.integrate_reference_z(r_anchored, slope, 6)
        @test za ≈ zb atol=1e-12
    end
end

@testset "lateral_offset_grid" begin
    @testset "straight line: pure perpendicular offset" begin
        x = [0.0, 1.0, 2.0, 3.0]   # straight along +x
        y = [0.0, 0.0, 0.0, 0.0]
        v = [-1.0, 0.0, 1.0]
        X, Y = OpenCRG.lateral_offset_grid(x, y, v)
        @test X ≈ repeat(x, 1, 3)                     # offsetting perpendicular to +x doesn't move x
        @test Y ≈ [-1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0] || Y ≈ -[-1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0]
        # (whichever sign convention `perp` uses; either is "correct" in isolation —
        # Task 17's cross-validation against the C library is the real arbiter.)
    end

    @testset "shape sanity: v=0 reproduces the reference line exactly" begin
        x = [0.0, 1.0, 2.5, 2.5]   # includes a kink, to exercise the miter-normal interior formula
        y = [0.0, 0.5, 1.0, 2.0]
        X, Y = OpenCRG.lateral_offset_grid(x, y, [0.0])
        @test X[:, 1] ≈ x
        @test Y[:, 1] ≈ y
    end

    @testset "interior miter join at nonzero v: a real, non-degenerate kink" begin
        # Neither test above actually exercises the interior miter-join formula's
        # correctness at a nonzero v. The "v=0" test above cannot: at v=0,
        # X[i,j] = x[i] + 0*offset_dir[i][1] == x[i] identically for ANY finite
        # offset_dir[i] whatsoever, however wrong -- v=0 multiplies the bisector
        # term away before it can affect the result (confirmed: offset_dir[2] for
        # that test's kink is actually (-0.372, 0.930), a nontrivial value that
        # the v=0 test never looks at). And the "straight line" test has no
        # interior kink at all (every segment shares the same direction, so
        # `chord` at each interior node equals `seg` itself and the whole
        # bisect-and-rescale computation degenerates to exactly each segment's
        # own normal) -- it can't distinguish a correct miter join from a wrong
        # one either. So this test suite, as given, never actually checks the
        # interior formula (`chord`/`n1`/`denom`/rescale) at a kink with v != 0.
        #
        # This closes that gap with a 3-node right-angle turn at NONZERO v.
        # Segment lengths are DELIBERATELY UNEQUAL (2, then 1) -- not just for
        # realism, but because equal segment lengths make the chord bisector
        # symmetric between the two flanking segments (n1 . n12[i] happens to
        # equal n1 . n12[i-1]), which would hide a bug that rescales against the
        # wrong (preceding, not following) neighbor's normal. Confirmed by
        # injecting that exact off-by-one bug (denom computed from n12[i-1]
        # instead of n12[i]) into a scratch copy of this function: with equal
        # segment lengths, buggy and correct offset_dir[2] were identical; with
        # the unequal lengths used below, they differ sharply ((-0.5, 1.0) buggy
        # vs (-1.0, 2.0) correct). A second scratch variant -- omitting the
        # rescale-by-denom entirely and using the raw unit bisector instead --
        # was also confirmed to diverge sharply from the values asserted below.
        #
        # Nodes: (0,0) -[length-2 segment, heading +x]-> (2,0)
        #             -[length-1 segment, heading +y]-> (2,1)   (a 90-degree left turn)
        # seg[1] = (1,0), seg[2] = (0,1); n12[1] = perp(seg[1]) = (0,1) and
        # n12[2] = perp(seg[2]) = (-1,0) are what boundary nodes 1 and 3 fall
        # back to directly. For the interior node 2:
        #   chord = normalize(x[3]-x[1], y[3]-y[1]) = normalize(2,1) = (2/sqrt(5), 1/sqrt(5))
        #   n1 = perp(chord) = (-1/sqrt(5), 2/sqrt(5))
        #   denom = n1 . n12[2] = (-1/sqrt(5))*(-1) + (2/sqrt(5))*0 = 1/sqrt(5)
        #   offset_dir[2] = n1 / denom = (-1, 2)        -- the sqrt(5) cancels exactly
        # Independently verified two ways: by hand (above) and in an isolated
        # Julia scratch session executing this exact chord/n1/denom formula,
        # both agreeing bit-for-bit (no floating-point residue at all, since
        # everything here cancels to exact integers).
        x = [0.0, 2.0, 2.0]
        y = [0.0, 0.0, 1.0]
        v = [-1.0, 1.0, 2.0]
        X, Y = OpenCRG.lateral_offset_grid(x, y, v)

        # offset_dir = [(0,1), (-1,2), (-1,0)]; X[i,j] = x[i] + v[j]*offset_dir[i][1], similarly for Y.
        # (These literals pin the CURRENT perp(dx,dy)=(-dy,dx) convention -- if a
        # later task flips that convention, these expected values would need
        # negating relative to x/y, same as the straight-line test above.)
        @test X ≈ [0.0 0.0 0.0; 3.0 1.0 0.0; 3.0 1.0 0.0]
        @test Y ≈ [-1.0 1.0 2.0; -2.0 2.0 4.0; 1.0 1.0 1.0]
    end

    @testset "exact U-turn: epsilon-guarded fallback, ported and cross-validated in Task 17" begin
        # Equal-length segments meeting at an exact 180-degree hairpin send the
        # miter rescale's denominator (and normalize2's chord length) to zero --
        # this isn't purely theoretical (it's the DEFAULT segment-length
        # configuration from integrate_reference_line's non-end-anchored branch).
        #
        # UPDATED during Task 17: this test originally pinned the UNGUARDED
        # behavior (`isnan(X[2,1])`) as a deliberate "documented gap, not yet
        # fixed" placeholder. Task 17's cross-validation loop added a synthetic
        # hairpin fixture (`synthetic_hairpin.crg`, binary/KDBI so its phi
        # channel encodes exact `Float64(pi)` with no ASCII-decimal truncation)
        # and confirmed the real C oracle does NOT produce NaN here: `crgEvaluv2xy`
        # returns the reference line's own point (to ~1e-16), for every v, at the
        # hairpin node -- i.e. `normalizeVector2`'s and the rescale step's
        # `1.0e-10` epsilon guards (crgEvaluv2xy.c ~lines 149/197) make the whole
        # bisector collapse gracefully to a ~zero vector instead of blowing up.
        # Those guards are now ported into `lateral_offset_grid` (see its
        # docstring), so this exact scenario -- where the chord skipped over the
        # hairpin node is bit-exactly `(0.0, 0.0)` -- now resolves to
        # `offset_dir[2] = (0.0, 0.0)` too (via `normalize2`'s `length < 1.0e-10`
        # branch returning its zero input unchanged, then the rescale's
        # `abs(denom) > 1.0e-10` check also failing on that zero `n1`), so the
        # offset point lands exactly ON the reference line, for any v.
        x = [0.0, 1.0, 0.0]
        y = [0.0, 0.0, 0.0]
        X, Y = OpenCRG.lateral_offset_grid(x, y, [1.0])
        @test X[2, 1] == 1.0   # == x[2], not merely ≈ -- the guarded path is exact here (chord is bit-exact zero)
        @test Y[2, 1] == 0.0   # == y[2]
    end
end

@testset "assemble_z_grid" begin
    r = OpenCRG.parse_road_crg(String[])
    z_grid = [0.0 0.1 0.2; 1.0 1.1 1.2]   # 2 rows x 3 v-columns
    z_ref = [10.0, 20.0]
    v = [-1.0, 0.0, 1.0]

    @testset "no banking: pure additive z_grid + z_ref" begin
        Z = OpenCRG.assemble_z_grid(z_grid, z_ref, nothing, r, v)
        @test Z ≈ [10.0 10.1 10.2; 21.0 21.1 21.2]
    end

    @testset "banking adds v * bank(u), clipped to [v_min, v_max]" begin
        banking = [0.05, -0.05]
        Z = OpenCRG.assemble_z_grid(z_grid, z_ref, banking, r, v)
        @test Z[1, :] ≈ [10.0 - 0.05, 10.1, 10.2 + 0.05]
        @test Z[2, :] ≈ [21.0 + 0.05, 21.1, 21.2 - 0.05]
    end

    @testset "banking clamp bounds come from the SAME v passed in -- self-referential no-op (see docstring)" begin
        # The plan's original version of this test asserted results clipped
        # against the OUTER-SCOPE `v = [-1.0, 0.0, 1.0]` defined above, as if
        # `assemble_z_grid` clipped its `wide_v` argument against some other,
        # independently-declared road-width axis. But the function has only
        # ONE `v` parameter, used BOTH to index z_grid's columns AND to
        # compute `vmin, vmax = first(v), last(v)` for the banking clamp --
        # so whatever `v` is actually passed in defines its own clamp
        # bounds; the outer `v` above is never seen by the function at all.
        #
        # Empirically confirmed by running the plan's test verbatim against
        # the plan's verbatim implementation before writing this fix -- it
        # FAILS both assertions:
        #   Z[1,1] ≈ z_grid[1,1] + z_ref[1] + 0.1*(-1.0)   # expected ≈ 9.9
        #   Evaluated: 9.5 ≈ 9.9  =>  false  (got 9.5)
        #   Z[1,3] ≈ z_grid[1,3] + z_ref[1] + 0.1*1.0       # expected ≈ 10.3
        #   Evaluated: 10.7 ≈ 10.299999999999999  =>  false  (got 10.7)
        # because vmin/vmax are actually computed from `wide_v` itself
        # (-5.0, 5.0), not from the outer-scope `v` -- so nothing is
        # actually clipped relative to `wide_v`'s own values.
        #
        # This isn't a one-off typo in the literals: it's a structural fact
        # for ANY sorted-ascending `v` (the only kind `CRGData.v` is ever
        # constructed as, via Task 9's `assemble_channels`'s `sortperm`) --
        # `first(v)`/`last(v)` are simply v's own min/max, so
        # `clamp(v[j], first(v), last(v)) == v[j]` for every `j`, always.
        # See `assemble_z_grid`'s docstring for why this differs from the C
        # reference (a continuous point-query evaluator, where the query v
        # and the channel's declared axis are genuinely two different
        # things) and why the clamp is kept anyway (harmless, defensive,
        # matches upstream intent, costs nothing).
        banking = [0.1, 0.1]
        wide_v = [-5.0, 0.0, 5.0]
        Z = OpenCRG.assemble_z_grid(z_grid, z_ref, banking, r, wide_v)
        @test Z[1, 1] ≈ z_grid[1,1] + z_ref[1] + 0.1 * (-5.0)   # NOT clipped: -5.0 IS wide_v's own min
        @test Z[1, 3] ≈ z_grid[1,3] + z_ref[1] + 0.1 * 5.0      # NOT clipped: 5.0 IS wide_v's own max
    end

    @testset "clamp arithmetic itself, via a synthetic non-monotonic v (never produced by the real pipeline)" begin
        # The previous testset shows NO sorted-ascending v can ever trigger
        # the clamp: v[1] <= v[j] <= v[end] holds by definition for any
        # ascending array, and `CRGData.v` is always sorted ascending. So
        # the clamp is mathematically inert on every real
        # `read_crg -> road_surface_grid` call path -- there is no
        # spec-compliant input that reaches it non-trivially.
        #
        # To still give the literal `vc = clamp(v[j], vmin, vmax)` line real
        # regression protection -- so that e.g. accidentally simplifying it
        # to `vc = v[j]` would be CAUGHT by the test suite, which the
        # no-op-on-realistic-v test above structurally cannot do -- this
        # feeds `assemble_z_grid` a deliberately non-monotonic `v`.
        # `assemble_z_grid` never asserts sortedness (it has no reason to,
        # since Task 9 already guarantees it upstream), so this is
        # syntactically accepted, even though `read_crg`/`assemble_channels`
        # would never actually hand it a `v` shaped like this one.
        #
        # vmin, vmax = first(v), last(v) = -1.0, 1.0 here, but the MIDDLE
        # entry v[2] = 5.0 lies outside that range, so it (and only it) gets
        # pulled down to vmax = 1.0; confirmed by direct execution, not just
        # by hand: Z[1,2] = z_grid[1,2] + z_ref[1] + 1.0*clamp(5.0,-1,1) =
        # 0.1 + 10.0 + 1.0*1.0 = 11.1 (the unclamped value would have been
        # 0.1 + 10.0 + 1.0*5.0 = 15.1 -- clearly different).
        banking = [1.0, 0.0]   # only row 1 is exercised below
        non_monotonic_v = [-1.0, 5.0, 1.0]
        Z = OpenCRG.assemble_z_grid(z_grid, z_ref, banking, r, non_monotonic_v)
        @test Z[1, 2] ≈ z_grid[1,2] + z_ref[1] + 1.0 * 1.0   # clamped down from 5.0 to vmax=1.0
        @test !isapprox(Z[1, 2], z_grid[1,2] + z_ref[1] + 1.0 * 5.0)   # sanity: distinguishable from the unclamped value
    end
end

using LibOpenCRG

@testset "apply_mods" begin
    @testset "no mods: identity (road_surface_grid unaffected)" begin
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_minimalist.crg"))
        d2 = OpenCRG.apply_mods(data)
        @test d2.phi == data.phi
        @test d2.refline.start_x == data.refline.start_x
    end

    @testset "REFLINE_OFFSET_*: rotate then translate, cross-validated against LibOpenCRG" begin
        path = joinpath(DATA, "handmade_curved_minimalist.crg")
        data = OpenCRG.read_crg(path)
        mods = OpenCRG.RoadCrgMods(refline_offset_phi=1.57, refline_offset_x=100.0, refline_offset_y=50.0)
        data_with_mods = OpenCRG.CRGData(data.comment, data.refline, data.format_code, data.opts, mods,
                                          data.mpro, data.phi, data.banking, data.slope, data.v, data.z)
        u, v, X, Y, Z = OpenCRG.road_surface_grid(data_with_mods)

        dsId = crgLoaderReadFile(path)
        # crgLoaderReadFile seeds dCrgModRefPointX/Y/Z/Phi defaults (via
        # crgOptionSetDefaultModifiers), and crgDataApplyTransformations checks
        # for ANY dCrgModRefPoint* entry BEFORE dCrgModRefLineOffset*, silently
        # short-circuiting the offset modifiers below if not cleared first --
        # this exact bug was found and fixed in Task 15's own test; see the
        # crgDataSetModifiersApply docstring in lib/LibOpenCRG/src/LibOpenCRG.jl.
        @test crgDataSetModifierRemoveAll(dsId) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetPhi, 1.57) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetX, 100.0) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetY, 50.0) != 0
        crgDataSetModifiersApply(dsId)
        cpId = crgContactPointCreate(dsId)
        for i in eachindex(u), j in eachindex(v)
            ref = crgEvaluv2xy(cpId, u[i], v[j])
            @test ref.status != 0
            @test X[i,j] ≈ ref.x atol=1e-6
            @test Y[i,j] ≈ ref.y atol=1e-6
        end
        crgContactPointDelete(cpId); crgDataSetRelease(dsId); crgMemRelease()
    end

    @testset "REFPOINT_* overrides REFLINE_OFFSET_*/REFLINE_ROTCENTER_* entirely, cross-validated" begin
        # apply_mods's docstring claims REFPOINT_* (if ANY such field is set)
        # takes over the rotate+translate step entirely, ignoring
        # REFLINE_OFFSET_*/REFLINE_ROTCENTER_* completely -- matching
        # `crgDataApplyTransformations` in crgMgr.c, where `applyXform` is set
        # by any of the `dCrgModRefPoint{X,Y,Z,Phi,U,UFrac,V,VFrac}` checks, and
        # the WHOLE REFLINE_OFFSET_*/ROTCENTER_* block only runs `if (!applyXform)`.
        # No test anywhere in this package actually exercised a REFPOINT_* field
        # before this one -- confirmed by reading crgDataApplyTransformations
        # directly (lib/LibOpenCRG/csrc/src/crgMgr.c, ~lines 804-922), then
        # closed here with a real cross-validation: refline_offset_*/rotcenter_*
        # are set to deliberately wild, easy-to-notice values (999.0, -999.0,
        # a 3.0 rad rotation, rotcenter far from the line) alongside a
        # refpoint_x/y/phi -- if has_refpoint's short-circuit were broken (e.g.
        # composing both instead of ignoring one), the result would be wildly
        # off from the C oracle, not subtly off.
        path = joinpath(DATA, "handmade_curved_minimalist.crg")
        data = OpenCRG.read_crg(path)
        mods = OpenCRG.RoadCrgMods(
            refpoint_x=10.0, refpoint_y=20.0, refpoint_phi=0.5,
            refline_offset_x=999.0, refline_offset_y=-999.0, refline_offset_phi=3.0,
            refline_rotcenter_x=-500.0, refline_rotcenter_y=500.0,
        )
        data_with_mods = OpenCRG.CRGData(data.comment, data.refline, data.format_code, data.opts, mods,
                                          data.mpro, data.phi, data.banking, data.slope, data.v, data.z)
        u, v, X, Y, Z = OpenCRG.road_surface_grid(data_with_mods)
        # refpoint_u/refpoint_u_fraction/refpoint_v/refpoint_v_fraction are all
        # unset, so the "from" point defaults to (u,v) = (start_u, 0) on BOTH
        # sides, matching crgDataApplyTransformations's own default (`uPos =
        # crgData->channelU.info.first; vPos = 0.0;`) -- so the reference line's
        # start point should land exactly on (refpoint_x, refpoint_y) with
        # heading refpoint_phi.
        d2 = OpenCRG.apply_mods(data_with_mods)
        @test d2.refline.start_x ≈ 10.0
        @test d2.refline.start_y ≈ 20.0
        @test d2.refline.start_phi ≈ 0.5

        dsId = crgLoaderReadFile(path)
        @test crgDataSetModifierRemoveAll(dsId) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefPointX, 10.0) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefPointY, 20.0) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefPointPhi, 0.5) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetX, 999.0) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetY, -999.0) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetPhi, 3.0) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineRotCenterX, -500.0) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineRotCenterY, 500.0) != 0
        crgDataSetModifiersApply(dsId)
        cpId = crgContactPointCreate(dsId)
        for i in eachindex(u), j in eachindex(v)
            ref = crgEvaluv2xy(cpId, u[i], v[j])
            @test ref.status != 0
            @test X[i,j] ≈ ref.x atol=1e-6
            @test Y[i,j] ≈ ref.y atol=1e-6
        end
        crgContactPointDelete(cpId); crgDataSetRelease(dsId); crgMemRelease()
    end

    @testset "SCALE_Z_GRID doubles elevation" begin
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        mods = OpenCRG.RoadCrgMods(scale_z_grid=2.0)
        data_with_mods = OpenCRG.CRGData(data.comment, data.refline, data.format_code, data.opts, mods,
                                          data.mpro, data.phi, data.banking, data.slope, data.v, data.z)
        d2 = OpenCRG.apply_mods(data_with_mods)
        # NOT plain `≈` here: this real fixture genuinely has 3 NaN samples
        # scattered in its z data (border columns at rows 12/14/16 -- confirmed
        # by direct inspection, `count(isnan, data.z) == 3`), and whole-array
        # `isapprox` is norm-based (`norm(d2.z - 2.0.*data.z) <= atol + ...`):
        # a NaN anywhere in the difference poisons the norm to NaN, which is
        # never `<=` anything, so the plan's original `d2.z ≈ 2.0 .* data.z`
        # fails even though every entry -- NaN included, since NaN*2.0 is still
        # NaN at the same position on both sides -- is exactly what it should
        # be. `isequal` is exact (not approximate) but that's fine: `.*=` here
        # is a single scalar multiply, bit-for-bit identical to `2.0 .* data.z`,
        # and `isequal`'s NaN-as-equal-to-itself semantics is exactly what's
        # needed to compare the missing-sample positions correctly.
        @test isequal(d2.z, 2.0 .* data.z)
    end

    @testset "SCALE_WIDTH preserves the v-ascending / z-column invariant assemble_z_grid relies on" begin
        # Task 14's review flagged that assemble_z_grid's banking clamp is only a
        # documented no-op because CRGData.v is always sorted ascending and lines
        # up 1:1 with CRGData.z's columns (Task 9's sortperm guarantee). apply_mods
        # is the only place in this task that touches `v` at all (via
        # SCALE_WIDTH), so this regression-tests that a (positive) width scale
        # can't desync that invariant: v must stay sorted ascending, and each
        # column of z must still correspond to the same (scaled) v it did before.
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        mods = OpenCRG.RoadCrgMods(scale_width=2.0)
        data_with_mods = OpenCRG.CRGData(data.comment, data.refline, data.format_code, data.opts, mods,
                                          data.mpro, data.phi, data.banking, data.slope, data.v, data.z)
        d2 = OpenCRG.apply_mods(data_with_mods)
        @test d2.v ≈ 2.0 .* data.v
        @test issorted(d2.v)
        # isequal, not ==: this fixture has real NaN samples (see the
        # SCALE_Z_GRID testset above) and plain `==` treats NaN != NaN even in
        # the exact-same-position, exact-same-untouched-value case here.
        @test isequal(d2.z, data.z)   # z itself is untouched -- only the v labels attached to its columns change
    end
end
