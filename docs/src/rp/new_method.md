# Contributing

## Adding a new rock physics equation

If you want to add a new rock physics equation, eg., a conductivity mechanism that got out recently and fits with your research or an anelasticity model you want to add, please open an issue/ pull request.

You can contribute by adding a new method following the steps listed below where we follow how a `conductivity` model is added but similar idea is used for `viscous`/`elastic`/`anelastic`:

 1. add a new `mutable struct` similar to other fomulations in `src/conductivity/types.jl`
 2. define `forward` dispatch in `src/conductivity/forward.jl`
 3. add the constants to be used in `forward` dispatch to `src/conductivity/cache.jl`. It is recommended to look at other `params` defined in the `cache.jl` file. Adding empirical constants in `params` keeps the codebase cleaner and organized.
 4. If required, add other functions in `src/conductivity/utils.jl`
 5. add a new distribution struct in `probabilistic/conductivity.jl` similar in definition as the `struct` defined in step 1.

It's always a good idea to add some docstring and test for the methods you add.
