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
end
