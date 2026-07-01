# lib/OpenCRG/test/runtests.jl
using Test
using OpenCRG

const DATA = joinpath(@__DIR__, "data")

@testset "OpenCRG.jl" begin
    include("test_header.jl")
    include("test_kddef.jl")
    include("test_payload.jl")
    include("test_read.jl")
    include("test_transform.jl")
    include("test_crossvalidate.jl")
end
