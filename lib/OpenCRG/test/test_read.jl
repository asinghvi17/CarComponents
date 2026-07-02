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
        @test !isempty(data.opts)   # this file has a real $ROAD_CRG_OPTS section -- guards against a swapped section-key typo
    end

    @testset "comment/opts/mods/mpro are wired to the right sections, not silently swapped" begin
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        @test data.mods isa OpenCRG.RoadCrgMods
        @test all(f -> getfield(data.mods, f) === nothing, fieldnames(OpenCRG.RoadCrgMods))
    end

    @testset "no channels declared is a clear error, not a DivideError" begin
        path = tempname()
        write(path, "\$CT\nempty file, no \$KD_DEFINITION at all\n\$\$\$\$\n")
        try
            @test_throws Exception OpenCRG.read_crg(path)
        finally
            rm(path)
        end
    end
end
