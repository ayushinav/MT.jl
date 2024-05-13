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

Looks like another change is coming up because when we want to do probabilistic inference over one field and not the other, the constraint that both the fields should be of the same type breaks because one will have `Dual` `eltype` while the other may have `Float64`.

What happens with precompilation is another detail we'll have to figure out soon...

No longer exporting `μ` because it is a common variable and may interact at places where `μ` is not defined properly. Better to encounter errors than running something with bugs.

**Ideation for abstract types**

Let's define the tree as :
We'll have a outermost abstract type `AbstractModel`, under which would be `AbstractGeophyModel`, make other subtypes as `AbstractGeoEMModel` but for now, let's add `MTModel` without `AbstractGeoEMModel`. 

Let's start by defining abstract type `AbstractResponse`, and then `AbstractGeophyResponse`, under which will be `MTResponse` for now.

For plotting, and/or any other purposes. we'll have `utils.jl`.

While we're at it, we should also specify more type information that the `MTmodel` field will acquire. We can also probably use `StaticArrays` inside occam codes, because all the vectors would be of the same size, similarly for MCMC models. Let's compare the performance for our current MT models. Also, know that if you're wanting to use `Float64`/ `Float32`, the autodiff packages require `Dual` types. Better to leave things at `AbstractArrays`.

Using parametric types under non parameteric abstract types is possible because the types do not inherit the fields of the abstract type unless done so explicitly. In short, types are for multiple dispatch and to make the codes informed about the types


Now, instead of having a 2D or a 3D type separately, we can simply do multiple dispatch on the `MTModel`s, eg., for 1d model, it would be `MTModel{Vector{T1}, Vector{T2}}`, for 2D, it would be `MTModel{Matrix{T1}, Matrix{T2}}`, for 3D, we would have `MTModel{Array{T1, 3}, Array{T2, 3}}`. SOMEWHAT similarly for `MTResponse`.