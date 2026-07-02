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
        # Backward pass (du=1): zb[5]=0.9, zb[4]=0.9-0.3=0.6, zb[3]=0.6-(-0.1)=0.7,
        # zb[2]=0.7-0.2=0.5, zb[1]=0.5-0.1=0.4 -> zb=[0.4, 0.5, 0.7, 0.6, 0.9].
        # Forward/blend recursion (chains off the already-blended z_ref[i], same
        # compounding structure as Task 11's (x,y) case, NOT plain linear
        # interpolation between 0.0 and 0.9):
        #   i=1, frac=1/4: z_ref[2] = 0.75*(0.0+0.1)  + 0.25*0.5 = 0.075  + 0.125 = 0.2
        #   i=2, frac=2/4: z_ref[3] = 0.5*(0.2+0.2)   + 0.5*0.7  = 0.2    + 0.35  = 0.55
        #   i=3, frac=3/4: z_ref[4] = 0.25*(0.55-0.1) + 0.75*0.6 = 0.1125 + 0.45  = 0.5625
        # Independently verified three ways: by hand, in an isolated Julia scratch
        # session, and via the deliberate bug-injection check described above.
        @test z_ref ≈ [0.0, 0.2, 0.55, 0.5625, 0.9] atol=1e-12
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
