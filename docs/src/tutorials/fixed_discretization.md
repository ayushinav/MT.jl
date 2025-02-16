## Fixed discretization
Geophysical models generally have fixed discretization. This is mostly because the different numerical schemes such as finite difference and finite element are computationally expensive and allocating a grid prior to solving the corresponding PDEs saves some computational resources. We provide the capability to do MCMC inference on such fixed grids.

Let's denote the model parameters, eg., conductivity, by `m`, and the layer thickness by `h`. Therefore, in a N-layer case, we will have 

```math
m = [m_1, m_2, m_3, ... , m_N] \\
h = [h_1, h_2, h_3, ... , h_{N_1}]
```

## Copy-Pasteable code
```julia
using MT
using Distributions
using Turing
using LinearAlgebra


m_test = MTModel(log10.([100., 10., 1000.]), [1e3, 1e3]);
f = 10 .^ range(-4, stop = 1, length = 25);
ω = vec(2π .* f);

r_obs = forward(m_test, ω);

err_phi = asin(0.01) * 180/π .* ones(length(ω));
err_appres = 0.02 * r_obs.ρₐ;
err_resp = MTResponse(err_appres, err_phi);

r_obs.ρₐ .= r_obs.ρₐ .+ err_appres;
r_obs.ϕ .= r_obs.ϕ .+ err_phi;

respD = MTResponseDistribution(normal_dist, normal_dist);

z = 10 .^collect(range(1, stop = 4, length = 100));
h = diff(z);

# fixed discretization
modelD = MTModelDistribution(
    Product(
    [Uniform(-1., 5.) for i in eachindex(z)]
    ),
    vec(h)
);


n_samples = 50;
mcache = mcmc_cache(modelD, respD, 50, NUTS());

mt_chain = stochastic_inverse(r_obs, err_resp, ω, mcache)
```

The obtained `mt_chain` contains the distributions that can be saved using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl).

```julia
using JLD2
JLD2.@save "file_path.jld2" mt_chain
```

**Note**: 
!!! note
    The returned chains will be sampled in the distribution specified by `modelD`. In the presented case, it will have values $\in [-1, 5]$ and we can get the values by `10. ^ value`.

The list of models can then be obtained from chains using
```
model_list = get_model_list(mt_chain, modelD)
```

We can then easily check the fit of the response curves
```
plt_resps = prepare_plot(r_obs, ω, alpha = 0.);
resp_models = forward(model_list[1], ω);

for i in 1:(length(model_list) > 50 ? 50 : length(model_list))
   forward!(resp_models, model_list[i], ω);
   prepare_plot!(resp_models, ω, alpha = 0.4); 
end

prepare_plot!(r_obs, ω, d_err = err_resp, markersize = 3, color = :orange);
plot_response(plt_resps)
```

The posterior distribution can then be obtained as:
```julia
pre_img = pre_image(m_dist, mt_chain);
kde_img = get_kde_image(pre_img..., false, xscale = :identity, yscale = :identity, yflip = true)
```

We can also obtain the mean and 1 std deviation bounds as:
```julia
mean_std_plt_lin = get_mean_std_image(pre_img..., yscale = :identity)
```