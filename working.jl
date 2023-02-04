using BenchmarkTools

include("fwd1d.jl")

T= 10 .^(range(-2,5,length= 57));
ω= 2π./T;
μ= 4π*1e-7;
# Z= zeros(size(T));

Z= fill(0. *im, length(T));
nω= length(T);

h= (2 .+randn(10)).* 5e2; # km
ρ= 10 .^(2 .+randn(11)); # Ωm

@mt_init nω;

@time for i in 1:1e6
    update_Z!(Z, ω, ρ, h);
end