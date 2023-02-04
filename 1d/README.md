1d MT forward modeling

Import all the functions and utilities


```julia
include("fwd1d.jl");
```

Set up the model space (here an 11-layer halfspace is initialized)


```julia
h= (2 .+randn(10)).* 5e2; # m
ρ= 10 .^(2 .+randn(11)); # Ωm
```

Set up the parameters for forward modeling


```julia
T= 10 .^(range(-2,5,length= 57));
ω= 2π./T;
μ= 4π*1e-7;

Z= fill(0. *im, length(T));
nω= length(T);
@mt_init nω; # initializes and allocates utilities for better performance
```

Test the speed!


```julia
using BenchmarkTools
```


```julia
@btime update_Z!(Z, ω, ρ, h);

```

      29.250 μs (0 allocations: 0 bytes)

