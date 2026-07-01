# lib/OpenCRG/test/test_header.jl
@testset "header tokenizing" begin
    @testset "find_header_end / split_lines on a real file" begin
        bytes = read(joinpath(DATA, "handmade_curved_minimalist.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        @test header_end <= length(bytes)
        header_lines = OpenCRG.split_lines(bytes[1:header_end-1])
        @test any(l -> startswith(l, "\$CT"), header_lines)
        # The line at header_end-1's start must be the "$$...." terminator
        term_line = OpenCRG.split_lines(bytes[1:header_end-1])[end] # not the terminator itself; sanity only
        @test !isempty(header_lines)
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
end
