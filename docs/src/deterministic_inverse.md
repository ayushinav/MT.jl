# Deterministic Inversion


## Brief introduction
Inverse problems in geophysics are notoriously ill-posed with non-unique solutions. MT inversion is no different. Below we demonstrate how we can perform different non-linear inverse schemes on a synthetic dataset.

## Demo
We start with defining models:

```@example inverse_demo
using MT, LinearAlgebra

ρ= log10.([100., 10., 400., 1000.])
h= [100., 1000., 10000.]
m= MTModel(ρ, h)

T= 10 .^(range(-1,5,length= 19))
ω= 2π./T

plot_model(m)
```

and getting data, with a 10% error floor.

```@example inverse_demo
resp= forward(m, ω)

err_resp = MTResponse(
    0.1 .* resp.ρₐ,
    asin(0.1) .+ zero(ω)
    )

plt= prepare_plot(resp, ω, err_resp, label= "observed")
plot_response(plt)
```

It's time to perform a few inversion schemes. All the inversion schemes can be called by the same function call, with just the difference of making an `alg_cache`, which mostly just depends on the library to be called for the non-linear inverse problem.

### Occam

Performing occam essentially boils down to making an `occam_cache`, which is done by making a call to [`Occam`](@ref).

```@example inverse_demo
W = diagm(inv.([err_resp.ρₐ..., err_resp.ϕ...])) .^ 2; # weight matrix

h_test= 10 .^range(0., 5., length= 50)
ρ_test= 2 .*ones(length(h_test)+1)

m_occam= MTModel(ρ_test, h_test);

alg_cache = Occam(; μgrid=[1e-2, 1e6])
inverse!(m_occam, resp, ω, alg_cache, W = W; max_iters= 50, verbose = true)
```

### Levenberg-Marquadt

While Occam is implemented in the package, we borrow a few from other packages. One of them is `NonlinearSolve.jl`, where we have `Levenberg-Marquadt` scheme. Again, everything boils down to creating the `alg_cache`. You can use [other solvers](https://docs.sciml.ai/NonlinearSolve/stable/solvers/nonlinear_least_squares_solvers/) from the package as well.

```@example inverse_demo
using NonlinearSolve

h_test= 10 .^range(0., 5., length= 50)
ρ_test= 2. .*ones(length(h_test)+1)

resp_trans_utils = (ρₐ=MT.log_tf, ϕ=MT.phi_scale_tf);

resp_lm = MTResponse(
    resp_trans_utils[:ρₐ].tf.(resp.ρₐ),
    resp_trans_utils[:ϕ].tf.(resp.ϕ),
)

err_resp_lm = MTResponse(
    resp_trans_utils[:ρₐ].dtf.(resp.ρₐ) .* err_resp.ρₐ,
    resp_trans_utils[:ϕ].dtf.(resp.ϕ) .* err_resp.ϕ,
)

m_lm= MTModel(ρ_test, h_test);

W_lm = diagm(inv.([err_resp_lm.ρₐ..., err_resp_lm.ϕ...])) .^ 2; # weight matrix

alg_cache =  NonlinearAlg(; alg = LevenbergMarquardt, μ = 1.0)
inverse!(m_lm, resp, ω, alg_cache; W =W_lm, max_iters= 100, verbose = true, response_trans_utils = (ρₐ=MT.log_tf, ϕ=MT.phi_scale_tf))
```

### LBFGS

Another popular algorithm is LBFGS, which we borrow from `Optimization.jl`. Again, create the `alg_cache` and it's good to go. `Optimization.jl` provides a [suite of solvers](https://docs.sciml.ai/Optimization/stable/#Overview-of-the-solver-packages-in-alphabetical-order), also by wrapping around a few others.

```@example inverse_demo
using Optimization, OptimizationOptimJL

h_test= 10 .^range(0., 5., length= 50)
ρ_test= m_occam.m .+ 1 .* randn(length(h_test) + 1)

resp_trans_utils = (ρₐ=MT.log_tf, ϕ=MT.phi_scale_tf);

resp_lbfgs = MTResponse(
    resp_trans_utils[:ρₐ].tf.(resp.ρₐ),
    resp_trans_utils[:ϕ].tf.(resp.ϕ),
)

err_resp_lbfgs = MTResponse(
    resp_trans_utils[:ρₐ].dtf.(resp.ρₐ) .* err_resp.ρₐ,
    resp_trans_utils[:ϕ].dtf.(resp.ϕ) .* err_resp.ϕ,
)

m_lbfgs= MTModel(ρ_test, h_test);

W_lbfgs = diagm(inv.([err_resp_lbfgs.ρₐ..., err_resp_lbfgs.ϕ...])) .^ 2; # weight matrix

alg_cache =  OptAlg(; alg = LBFGS, μ = 1.0)
inverse!(m_lbfgs, resp, ω, alg_cache, W = W_lbfgs; max_iters= 50, verbose = true, response_trans_utils = (ρₐ=MT.log_tf, ϕ=MT.phi_scale_tf))
```

## Fits

So how well do we fit the data? Note that, in no way we compare the different inversion schemes here. A lot of these schemes depend heavily on the initial model and our choice might be sub-optimal.

```@example inverse_demo
resp_occam= forward(m_occam, ω);
resp_lm= forward(m_lm, ω);
resp_lbfgs= forward(m_lbfgs, ω);

plt= prepare_plot(resp, ω, err_resp, label= "true", plt_type = :scatter)
prepare_plot!(plt, resp_occam, ω, label= "occam", plt_type = :plot, linewidth= 1)
prepare_plot!(plt, resp_lm, ω, label= "Levenberg-Marquadt", plt_type = :plot, linewidth= 1)
prepare_plot!(plt, resp_lbfgs, ω, label= "LBFGS", plt_type = :plot, linewidth= 1)
plot_response(plt,legend=:bottomright)
```

And a look at different models

```@example inverse_demo
plt = plot_model(m, label = "true", linewidth = 3, color = "black")
plot_model!(plt, m_occam, label = "occam", linewidth = 2, color = "blue")
plot_model!(plt, m_lm, label = "Levenberg-Marquadt", linewidth = 2)
plot_model!(plt, m_lbfgs, label = "LBFGS", linewidth = 2)
```