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

    # LDFI closes a whole-plan-final-review coverage gap: despite being one of
    # the 4 "must parse" encodings, grep confirmed :LDFI was never actually
    # decoded by any test before this fix -- it only ever appeared in the
    # ASCII_FIELD_WIDTH/FIELDS_PER_RECORD Dict literals above and docstring
    # prose. Unlike LRFI/KRBI's real vendored fixtures, no upstream example
    # file uses LDFI, so (mirroring the existing KDBI precedent immediately
    # below in the binary section, "no upstream fixture covers this") these
    # are synthetic, built directly against decode_ascii_payload rather than
    # a full read_crg round-trip -- cheaper, and matches this file's existing
    # unit-test granularity for the other format lacking a real fixture.
    # Fields are built with `lpad(..., 20)` rather than hand-typed fixed-width
    # literals -- LDFI's field width (20) is exactly double LRFI's (10), and
    # hand-counting spaces is exactly the kind of silent-misalignment mistake
    # the LRFI tests above warn about; lpad guarantees the width by
    # construction instead of by careful counting.
    @testset "decode_ascii_payload: exact multiple of per_record, no line wrap (4 channels at LDFI)" begin
        vals = [-1.0, 0.0, 1.0, 2.5]
        line = join(lpad(string(v), 20) for v in vals)
        @test length(line) == 80   # 4 fields x 20 chars, sanity-checking the synthesized fixture itself
        m = OpenCRG.decode_ascii_payload([line], :LDFI, 4)
        @test size(m) == (1, 4)
        @test m[1, :] == vals
    end

    @testset "decode_ascii_payload: row-wrapping, 6 channels at LDFI (4/record)" begin
        # 6 channels at 4 fields/record wraps onto 2 lines (4 then 2) -- the
        # LDFI-specific analogue of the LRFI row-wrapping test above, using a
        # DIFFERENT per-record count (4, not 8) so this can't silently pass by
        # accidentally reusing LRFI's per-record constant.
        row = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
        line1 = join(lpad(string(v), 20) for v in row[1:4])
        line2 = join(lpad(string(v), 20) for v in row[5:6])
        m = OpenCRG.decode_ascii_payload([line1, line2], :LDFI, 6)
        @test size(m) == (1, 6)
        @test m[1, :] == row
    end

    @testset "decode_ascii_payload: **unused** literal at LDFI's actual 20-char field width" begin
        # Task 7's docstring (decode_ascii_field, above) flags this as subtle:
        # the bare 10-char literal "**unused**" structurally can never
        # string-equal a raw 20-char LDFI field, so it must fall through to
        # the numeric-parse path and decode to NaN, not 0.0 -- unlike LRFI,
        # where the same literal (exactly 10 chars) matches directly and
        # decodes to 0.0. The existing decode_ascii_field unit test above
        # (`"**unused**" * " "^10`) already pins this at the field-decoding
        # level in isolation; this test additionally exercises it through
        # decode_ascii_payload's REAL fixed-width slicing at nchannels=4,
        # confirming the 20-char slice taken from a realistic multi-field
        # line lands on the right offsets -- i.e. that neighboring numeric
        # fields decode correctly (proving the NaN field didn't desync the
        # per-field width bookkeeping), not just that the literal alone
        # decodes to NaN when handed to decode_ascii_field directly.
        unused_field = "**unused**" * " "^10
        @test length(unused_field) == 20
        line = lpad("3.5", 20) * unused_field * lpad("-2.25", 20) * lpad("100.0", 20)
        m = OpenCRG.decode_ascii_payload([line], :LDFI, 4)
        @test m[1, 1] == 3.5
        @test isnan(m[1, 2])          # **unused**, padded to LDFI width -- NaN, not 0.0
        @test m[1, 3] == -2.25        # confirms field 2's NaN didn't shift field 3's offset
        @test m[1, 4] == 100.0        # ...nor field 4's
    end
end

@testset "binary payload decoding" begin
    @testset "synthetic 2-row x 3-channel KRBI blob catches row/column transposition" begin
        # Row 0 = [1.0, 2.0, 3.0], Row 1 = [4.0, 5.0, 6.0], packed tightly,
        # big-endian Float32, no padding between rows (unlike ASCII).
        vals = Float32[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        bytes = UInt8[]
        for v in vals
            append!(bytes, reverse(reinterpret(UInt8, [v])))  # host is little-endian -> reverse for big-endian
        end
        m = OpenCRG.decode_binary_payload(bytes, :KRBI, 3)
        @test size(m) == (2, 3)
        @test m[1, :] == [1.0, 2.0, 3.0]
        @test m[2, :] == [4.0, 5.0, 6.0]
    end

    @testset "synthetic 2-row x 3-channel KDBI blob (Float64, no upstream fixture covers this)" begin
        vals = Float64[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        bytes = UInt8[]
        for v in vals
            append!(bytes, reverse(reinterpret(UInt8, [v])))
        end
        m = OpenCRG.decode_binary_payload(bytes, :KDBI, 3)
        @test size(m) == (2, 3)
        @test m[1, :] == [1.0, 2.0, 3.0]
        @test m[2, :] == [4.0, 5.0, 6.0]
    end

    @testset "truncated row is an error, not a silent drop" begin
        row = Float32.(1:25)
        bytes = UInt8[]
        for v in row
            append!(bytes, reverse(reinterpret(UInt8, [v])))
        end
        partial_second_row = bytes[1:90]   # 90 leftover bytes: not plausible padding (>= 80)
        @test_throws Exception OpenCRG.decode_binary_payload(vcat(bytes, partial_second_row), :KRBI, 25)
    end

    @testset "real binary file: shape and no per-row 80-byte alignment" begin
        bytes = read(joinpath(DATA, "belgian_block.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        payload = bytes[header_end:end]
        m = OpenCRG.decode_binary_payload(payload, :KRBI, 342)
        @test size(m) == (1001, 342)   # nu is DERIVED: 1369440 bytes ÷ 1368 bytes/row = 1001 (72 trailing padding bytes ignored)
        @test isnan(m[1, 1])                    # row 0 phi placeholder, per spec/research
        @test m[2, 1] ≈ 2.6527974605560303       # bit-exact match to REFERENCE_LINE_START_PHI in this file's header
        @test all(isnan, m[1, 2:22])             # first cross-section's missing left-border samples (verified byte-for-byte: channels 2-22 are NaN)
        @test isfinite(m[1, 23])                 # first real elevation sample in row 0 (verified: channel 23 == 2.1214599609375)
    end
end

@testset "assemble_channels" begin
    r = OpenCRG.parse_road_crg(["REFERENCE_LINE_START_PHI = 0.25"])
    channels = [
        OpenCRG.ChannelDef(:phi, nothing, nothing),
        OpenCRG.ChannelDef(:banking, nothing, nothing),
        OpenCRG.ChannelDef(:long_section, 1.0, nothing),   # declared out of ascending order on purpose
        OpenCRG.ChannelDef(:long_section, -1.0, nothing),
    ]
    raw = [  # 2 rows x 4 channels, in DECLARATION order (phi, banking, v=1.0, v=-1.0)
        99.0  0.1  10.0  20.0
        0.5   0.2  11.0  21.0
    ]
    phi, banking, slope, v, z = OpenCRG.assemble_channels(raw, channels, r)
    @test phi[1] == 0.25            # row-0 placeholder overwritten with REFERENCE_LINE_START_PHI...
    @test phi[2] == 0.5             # ...but row 1 keeps its real stored value
    @test banking == [0.1, 0.2]
    @test slope === nothing
    @test v == [-1.0, 1.0]          # sorted ascending, regardless of declaration order
    @test z == [20.0 10.0; 21.0 11.0]   # columns reordered to match the sorted v
end

@testset "assemble_channels: no long_section channels is an error" begin
    r = OpenCRG.parse_road_crg(String[])
    channels = [OpenCRG.ChannelDef(:phi, nothing, nothing)]
    raw = reshape([0.0, 0.5], 2, 1)
    @test_throws Exception OpenCRG.assemble_channels(raw, channels, r)
end

@testset "assemble_channels: no phi channel declared -- NaN propagates, doesn't error" begin
    r = OpenCRG.parse_road_crg(["REFERENCE_LINE_START_PHI = 0.25"])
    channels = [OpenCRG.ChannelDef(:long_section, 0.0, nothing)]
    raw = reshape([10.0, 11.0], 2, 1)
    phi, banking, slope, v, z = OpenCRG.assemble_channels(raw, channels, r)
    @test phi[1] == 0.25      # row-1 placeholder still gets overwritten...
    @test isnan(phi[2])       # ...but there's no real data to recover the rest from
end
