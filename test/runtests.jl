using MT
using Test

@testset "MT.jl" begin
    h= (2 .+randn(10)).* 5e2; # m
    ρ= 10 .^(2 .+randn(11)); # Ωm
    T= 10 .^(range(-2,5,length= 57));
    ω= 2π./T;

    Z= fill(0. *im, length(T));
    nω= length(T);
    @mt_init nω; #
    update_Z!(Z, ω, ρ, h);

    # Write necessary tests
end
