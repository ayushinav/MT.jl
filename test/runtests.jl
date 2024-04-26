using MT
using LinearAlgebra
using Test
using BenchmarkTools

@testset "MT.jl" begin
    include("forward_test.jl")
    include("inverse_test.jl")
end