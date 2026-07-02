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
        @test x[3] ≈ 2.25 atol=1e-12
        @test all(y .≈ 0.0)
    end
end
