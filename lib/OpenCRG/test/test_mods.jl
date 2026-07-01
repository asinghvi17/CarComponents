# lib/OpenCRG/test/test_mods.jl
@testset "\$ROAD_CRG_MODS parsing" begin
    lines = [
        "SCALE_Z_GRID = 2.0",
        "REFLINE_OFFSET_PHI = 1.57",
        "REFLINE_OFFSET_X = 100.0",
        "REFPOINT_PHI = 0.0",   # presence alone should flip has_refpoint on, even though it's 0.0
    ]
    m = OpenCRG.parse_road_crg_mods(lines)
    @test m.scale_z_grid == 2.0
    @test m.scale_slope === nothing
    @test m.refline_offset_phi == 1.57
    @test m.refline_offset_x == 100.0
    @test m.refpoint_phi == 0.0

    empty_mods = OpenCRG.parse_road_crg_mods(String[])
    @test empty_mods.scale_z_grid === nothing
    @test empty_mods.refpoint_phi === nothing

    @testset "every RoadCrgMods field round-trips through its own header key" begin
        # NOTE: unrolled into one @test per field (rather than looping over
        # field_to_key and calling `getfield(m, field) == sentinel`) because the
        # looped form's failure output doesn't name the offending field: Julia's
        # @test only decomposes the top-level `==` call, printing e.g.
        # "Evaluated: nothing == 12345.6789" without revealing what `field` was.
        # Writing out `m.<field> == sentinel` per field makes the failing
        # expression itself name the field, so a typo'd key is immediately
        # attributable to one line.
        field_to_key = Dict(
            :scale_z_grid => "SCALE_Z_GRID", :scale_slope => "SCALE_SLOPE", :scale_banking => "SCALE_BANKING",
            :scale_length => "SCALE_LENGTH", :scale_width => "SCALE_WIDTH", :scale_curvature => "SCALE_CURVATURE",
            :grid_nan_offset => "GRID_NAN_OFFSET",
            :refpoint_u => "REFPOINT_U", :refpoint_u_fraction => "REFPOINT_U_FRACTION", :refpoint_u_offset => "REFPOINT_U_OFFSET",
            :refpoint_v => "REFPOINT_V", :refpoint_v_fraction => "REFPOINT_V_FRACTION", :refpoint_v_offset => "REFPOINT_V_OFFSET",
            :refpoint_x => "REFPOINT_X", :refpoint_y => "REFPOINT_Y", :refpoint_z => "REFPOINT_Z", :refpoint_phi => "REFPOINT_PHI",
            :refline_rotcenter_x => "REFLINE_ROTCENTER_X", :refline_rotcenter_y => "REFLINE_ROTCENTER_Y",
            :refline_offset_x => "REFLINE_OFFSET_X", :refline_offset_y => "REFLINE_OFFSET_Y",
            :refline_offset_z => "REFLINE_OFFSET_Z", :refline_offset_phi => "REFLINE_OFFSET_PHI",
        )
        sentinel = 12345.6789
        lines = ["$(key) = $(sentinel)" for key in values(field_to_key)]
        m = OpenCRG.parse_road_crg_mods(lines)
        @test m.scale_z_grid == sentinel
        @test m.scale_slope == sentinel
        @test m.scale_banking == sentinel
        @test m.scale_length == sentinel
        @test m.scale_width == sentinel
        @test m.scale_curvature == sentinel
        @test m.grid_nan_offset == sentinel
        @test m.refpoint_u == sentinel
        @test m.refpoint_u_fraction == sentinel
        @test m.refpoint_u_offset == sentinel
        @test m.refpoint_v == sentinel
        @test m.refpoint_v_fraction == sentinel
        @test m.refpoint_v_offset == sentinel
        @test m.refpoint_x == sentinel
        @test m.refpoint_y == sentinel
        @test m.refpoint_z == sentinel
        @test m.refpoint_phi == sentinel
        @test m.refline_rotcenter_x == sentinel
        @test m.refline_rotcenter_y == sentinel
        @test m.refline_offset_x == sentinel
        @test m.refline_offset_y == sentinel
        @test m.refline_offset_z == sentinel
        @test m.refline_offset_phi == sentinel

        # grid_nan_mode is Int-typed, not Float64 -- test separately with an integral value
        m2 = OpenCRG.parse_road_crg_mods(["GRID_NAN_MODE = 3"])
        @test m2.grid_nan_mode == 3
    end
end
