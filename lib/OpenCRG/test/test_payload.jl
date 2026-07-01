# lib/OpenCRG/test/test_payload.jl
@testset "ASCII payload decoding" begin
    @testset "decode_ascii_field" begin
        @test OpenCRG.decode_ascii_field("**unused**") == 0.0   # exact literal only
        @test isnan(OpenCRG.decode_ascii_field("*missing*"))
        @test OpenCRG.decode_ascii_field(" 0.0111111") == 0.0111111
        @test OpenCRG.decode_ascii_field("-1.500000") == -1.5
        @test OpenCRG.decode_ascii_field("1.0D+02") == 100.0   # Fortran D-exponent normalization
        @test isnan(OpenCRG.decode_ascii_field("**unused**" * " "^10))   # 20-char LDFI-width field: never matches the 10-char literal
    end

    @testset "decode_ascii_payload: exact multiple of per_record, no line wrap (8 channels at LRFI)" begin
        # Each LRFI field is a fixed 10-char width, right-justified -- NOT merely
        # space-separated 6-decimal numbers (a naive single-space join is only 9 chars wide
        # for non-negative values like "0.000000", so it silently misaligns every field after
        # the first; verified directly by round-tripping through decode_ascii_payload).
        lines = [" -1.000000  0.000000  1.000000  2.000000  3.000000  4.000000  5.000000  6.000000"]
        m = OpenCRG.decode_ascii_payload(lines, :LRFI, 8)
        @test size(m) == (1, 8)
        @test m[1, :] == [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    end

    @testset "decode_ascii_payload: non-whole-row remainder is an error, not a silent truncation" begin
        lines = [
            "**unused** 0.0000000**unused** 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000",
            " 0.0000000 0.0000000",
            " 0.0000000 0.0000000",   # a stray 3rd line: 3 lines don't divide evenly by 2 lines/row
        ]
        @test_throws Exception OpenCRG.decode_ascii_payload(lines, :LRFI, 10)
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
        # Channel 3 is "slope" (not "banking" -- the file's channel order is phi, banking,
        # slope, then 7 long-sections). At row 0 it is genuinely the exact literal
        # "**unused**" in the real file (confirmed by reading the raw bytes), so it decodes
        # to 0.0, not NaN.
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
