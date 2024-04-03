We move from 

```julia
mutable struct model{T <: Union{AbstractVector{Float32}, AbstractVector{Float64}}}
    m::T
    h::T
end
```
to

```julia
mutable struct model{T}
    m::T
    h::T
end
```

for both `model`and `response`. We probably started. For quite some time, I was confused if defining the types (`Vector{Float64}`) would shrink down the precompilation time. The difference is not much. This code is more abstract, allows for using different libraries, eg., the earlier version errored when using `Turing.jl` because it requires `Dual` numbers to get passed whereas we had a function only for `Float`s. This also opens the doors for precision. If all the operation are on `Float32`, after precompilation, we will have a specialized function for `Float32`. The constant `μ` will have type promotion and it's very rare that anything below `Float32` precision will get passed.

No longer exporting `μ` because it is a common variable and may interact at places where `μ` is not defined properly. Better to encounter errors than running something with bugs