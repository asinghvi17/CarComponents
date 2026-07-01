# lib/OpenCRG/test/test_kddef.jl
@testset "\$KD_DEFINITION parsing" begin
    @testset "position-form (at v = ...), U: lines ignored" begin
        lines = [
            "#:LRFI",
            "U:reference line u,m,0,1.0",
            "D:reference line phi,rad",
            "D:reference line banking,m/m",
            "D:reference line slope,m/m",
            "D:long section at v = -1.500,m",
            "D:long section at v =  0.000,m",
            "D:long section at v =  1.500,m",
        ]
        format_code, channels = OpenCRG.parse_kd_definition(lines)
        @test format_code == :LRFI
        @test length(channels) == 6   # the U: line contributes nothing
        @test [c.kind for c in channels] == [:phi, :banking, :slope, :long_section, :long_section, :long_section]
        @test [c.v for c in channels if c.kind == :long_section] == [-1.5, 0.0, 1.5]
    end

    @testset "index-form (bare N), default format code" begin
        lines = ["D:long section 1,m", "D:long section 2,m", "D:long section 3,m"]
        format_code, channels = OpenCRG.parse_kd_definition(lines)
        @test format_code == :KRBI   # default when no #: line is present
        @test [c.index for c in channels] == [1, 2, 3]
    end

    @testset "mixed position/index form is an error" begin
        lines = ["D:long section at v = 0.0,m", "D:long section 2,m"]
        @test_throws Exception OpenCRG.parse_kd_definition(lines)
    end

    @testset "v_axis: index-form uniform spacing" begin
        r = OpenCRG.parse_road_crg(["LONG_SECTION_V_RIGHT = -1.0", "LONG_SECTION_V_LEFT = 1.0", "LONG_SECTION_V_INCREMENT = 1.0"])
        channels = [OpenCRG.ChannelDef(:long_section, nothing, i) for i in 1:3]
        @test OpenCRG.v_axis(channels, r) == [-1.0, 0.0, 1.0]
    end

    @testset "v_axis: index-form derives increment from v_right/v_left when not explicitly given" begin
        r = OpenCRG.parse_road_crg(["LONG_SECTION_V_RIGHT = -1.0", "LONG_SECTION_V_LEFT = 1.0"])
        channels = [OpenCRG.ChannelDef(:long_section, nothing, i) for i in 1:5]
        @test OpenCRG.v_axis(channels, r) == [-1.0, -0.5, 0.0, 0.5, 1.0]
    end

    @testset "v_axis: position-form reads channel v directly, including non-uniform spacing" begin
        r = OpenCRG.parse_road_crg(String[])
        channels = [OpenCRG.ChannelDef(:long_section, v, nothing) for v in [-1.0, 0.0, 0.3]]
        @test OpenCRG.v_axis(channels, r) == [-1.0, 0.0, 0.3]
    end

    @testset "real file: handmade_curved_banked_sloped.crg" begin
        bytes = read(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        sections = OpenCRG.group_sections(OpenCRG.split_lines(bytes[1:header_end-1]))
        format_code, channels = OpenCRG.parse_kd_definition(sections["KD_DEFINITION"])
        @test format_code == :LRFI
        @test length(channels) == 10   # phi, banking, slope, 7 long sections
    end
end
