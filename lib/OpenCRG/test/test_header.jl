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
end
