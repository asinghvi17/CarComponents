# lib/OpenCRG/test/test_header.jl
@testset "header tokenizing" begin
    @testset "find_header_end / split_lines on a real file" begin
        bytes = read(joinpath(DATA, "handmade_curved_minimalist.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        @test header_end <= length(bytes)
        header_lines = OpenCRG.split_lines(bytes[1:header_end-1])
        @test any(l -> startswith(l, "\$CT"), header_lines)
        @test !isempty(header_lines)
        @test all(c -> c == '$', header_lines[end])   # the last header line is the "$$..." terminator
    end

    @testset "group_sections strips comments and groups by keyword" begin
        lines = [
            "\$ROAD_CRG",
            "* this whole line is a comment",
            "REFERENCE_LINE_INCREMENT = 1.0   ! inline comment here",
            "\$KD_DEFINITION",
            "#:LRFI",
        ]
        sections = OpenCRG.group_sections(lines)
        @test sections["ROAD_CRG"] == ["REFERENCE_LINE_INCREMENT = 1.0"]
        @test sections["KD_DEFINITION"] == ["#:LRFI"]
    end

    @testset "binary file: header/payload boundary is byte-exact" begin
        bytes = read(joinpath(DATA, "belgian_block.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        payload = bytes[header_end:end]
        # 1001 rows x 342 channels x 4 bytes, rounded up to a multiple of 80
        expected_padded = cld(1001 * 342 * 4, 80) * 80
        @test length(payload) == expected_padded
    end

    @testset "group_sections on a real file never produces an all-dollar section key" begin
        bytes = read(joinpath(DATA, "belgian_block.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        sections = OpenCRG.group_sections(OpenCRG.split_lines(bytes[1:header_end-1]))
        @test all(k -> !all(==('$'), k), keys(sections))
    end

    @testset "parse_road_crg on a real file" begin
        bytes = read(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        sections = OpenCRG.group_sections(OpenCRG.split_lines(bytes[1:header_end-1]))
        r = OpenCRG.parse_road_crg(sections["ROAD_CRG"])
        @test r.start_u == 0.0
        @test r.end_u == 22.0
        @test r.increment == 1.0
        @test r.start_y == 0.0
        @test r.start_phi == 0.0
        @test r.end_x === nothing   # this file has no explicit end position
        @test r.end_phi == 0.0   # present-but-zero must stay a real Float64, not become `nothing`
        @test r.v_right == -1.5
        @test r.v_left == 1.5
        @test r.start_slope == 0.0   # default when REFERENCE_LINE_START_S is absent
        @test r.start_banking == 0.0
    end

    @testset "parse_keyvalues / parse_keyvalue_strings" begin
        @test OpenCRG.parse_keyvalues(["FOO = 1.5", "BAR=2"]) == Dict("FOO"=>1.5, "BAR"=>2.0)
        @test OpenCRG.parse_keyvalue_strings(["PROJ_NM = UTM", "PROJ_ZONE = 32"]) ==
            Dict("PROJ_NM"=>"UTM", "PROJ_ZONE"=>"32")
    end

    @testset "\$ROAD_CRG_OPTS / \$ROAD_CRG_MPRO parse as generic dicts (not applied)" begin
        opts = OpenCRG.parse_keyvalues(["BORDER_MODE_U = 2", "BORDER_MODE_V = 0"])
        @test opts["BORDER_MODE_U"] == 2.0

        mpro = OpenCRG.parse_keyvalue_strings(["PROJ_NM = UTM", "GELL_A = 6378137.0"])
        @test mpro["PROJ_NM"] == "UTM"
        @test mpro["GELL_A"] == "6378137.0"   # kept as a string; no geodesy math is implemented
    end
end
