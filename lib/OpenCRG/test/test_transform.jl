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
