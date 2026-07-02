# lib/OpenCRG/test/test_read.jl
@testset "read_crg" begin
    @testset "small ASCII file end-to-end" begin
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_minimalist.crg"))
        @test data isa OpenCRG.CRGData
        @test data.format_code == :LRFI
        @test length(data.phi) == 23      # derived from the payload's actual row count — this file has no REFERENCE_LINE_END_U at all
        @test size(data.z, 1) == 23
        @test size(data.z, 2) == length(data.v)
    end

    @testset "banked/sloped ASCII file has banking and slope channels" begin
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        @test data.banking !== nothing
        @test data.slope !== nothing
        @test length(data.v) == 7
    end

    @testset "real binary file end-to-end" begin
        data = OpenCRG.read_crg(joinpath(DATA, "belgian_block.crg"))
        @test data.format_code == :KRBI
        @test size(data.z) == (1001, 341)   # 342 channels total minus 1 phi channel
    end
end
