# lib/OpenCRG/test/test_payload.jl
@testset "ASCII payload decoding" begin
    @testset "decode_ascii_field" begin
        @test OpenCRG.decode_ascii_field("**unused**") == 0.0   # exact literal only
        @test isnan(OpenCRG.decode_ascii_field("*missing*"))
        # decode_ascii_field's own docstring documents its precondition: the input is
        # "already stripped of surrounding whitespace by the caller" (decode_ascii_payload
        # always calls it as decode_ascii_field(strip(field))). Calling it directly with an
        # unstripped, space-padded field -- as fixed-width LRFI fields commonly are -- would
        # violate that precondition and spuriously return NaN (the space isn't in the accepted
        # character set), so this test strips first, matching real call-site usage.
        @test OpenCRG.decode_ascii_field(strip(" 0.0111111")) == 0.0111111
        @test OpenCRG.decode_ascii_field("-1.500000") == -1.5
    end

    @testset "decode_ascii_payload: row-wrapping, 10 channels at LRFI (8/record)" begin
        # Row 0 from handmade_curved_banked_sloped.crg, 10 channels, wraps onto 2 lines.
        lines = [
            "**unused** 0.0000000**unused** 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000",
            " 0.0000000 0.0000000",
        ]
        m = OpenCRG.decode_ascii_payload(lines, :LRFI, 10)
        @test size(m) == (1, 10)
        @test m[1,1] == 0.0          # **unused** -> 0.0
        @test m[1,2] == 0.0
        # Channel 3 (slope) at row 0 is, byte-for-byte in the real fixture, the exact
        # literal "**unused**" (verified directly against handmade_curved_banked_sloped.crg
        # rather than assumed) -- NOT "*missing*". decode_ascii_field only maps the exact
        # 10-char literal "**unused**" to 0.0, so this decodes to 0.0, not NaN.
        @test m[1,3] == 0.0
    end

    @testset "real file end-to-end shape" begin
        bytes = read(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        payload_lines = OpenCRG.split_lines(bytes[header_end:end])
        m = OpenCRG.decode_ascii_payload(payload_lines, :LRFI, 10)
        @test size(m) == (23, 10)    # nu is DERIVED from payload size: 46 lines / 2 lines-per-row = 23
        @test all(isfinite, m[2:end, 1])   # phi is defined for every row except row 1 (index 1 in Julia)
    end
end
