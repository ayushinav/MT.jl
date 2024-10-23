<!-- # Probabilistic inversion

We use [Turing.jl](https://github.com/TuringLang/Turing.jl) to perform stochastic inversion. 

## Constructing distributions to sample from
### Model distribution (*a priori* information)
Before beginning to talk about how to construct the *a priori* distribution, it is important to understand that the model here refers to the discretization space as well as the values of physical properties. For eg, for electrical methods, the model space will consist of the electrical resistivity as well as the grid sizes represented by those. This is primarily useful in 1D, but can have its own consequences.

For most applications, however, we fix the node points. We follow a 1D MT example to show the framework. For now, we begin by choosing a very broad prior for all the `n` layers. A model distribution can be constructed using `MTModelDistribution(...)`. The first parameter denotes the prior for electrical conductivities, while the second is for the layer thicknesses `h`. If you do not want to infer on `h`, just pass it as a simple vector, as is also demonstrated in the following case:

```julia
n = 50 # number of layers
h = collect(range(start = 0, stop = 1e3, length = n)) # making the discretization

modelD = MTModelDistribution(
    Product(
    [Uniform(-1, 5) for i in 1:n]
    ),
    vec(h)
);
```

Some of the things to understand here are :
* `Product(...)` is a `Distributions.jl` function that joins a series of multivariate samples, here `Uniform(-1,5)`. Thus, `Product([Uniform(-1, 5)])` makes a sampler that will output an `n`-length vector. Here, we use a uniform prior, but one can choose any prior. The only thing to keep in mind is that, the prior can be sampled. Another eg., would be `MultivariateNormal(μ, Σ)`.
* Notice that we're sampling on the log-scale, but haven't put that information out here yet. This transformation will go into the function we'll call for inverison.
* To also sample `h`, pass another multivariate distribution that outputs `n-1` length vector. Again using a uniform, independent prior, we'll have:
```julia
modelD = MTModelDistribution(
    Product(
    [Uniform(-1, 5) for i in 1:n]
    ),
    Product(
    [Uniform(50, 100) for i in 1:n-1]
    )
);
```
This has its own consequences on the results and one needs to be mindful of that while interpreting the results.

### Response distribution (*likelihood*)
A likelihood is determined by an observed response `r_obs` we want to fit, and the errors associated with it `err_resp`. One of the popular ways in which likelihood can be formed is using the gaussian distribution, centered around `r_obs` with variance in `err_resp`. To make things consistent, we pass `r_obs` and `err_resp` into the final function. To make a likelihood, we just need to pass a function that can take in a response parameter and the associated error and produce a distribution. We already provide a function `norm_dist` that takes in a vector for the response and another vector/matrix for the covariance matrix of the errors. The response distribution is then constructed by:
```julia
respD = MTResponseDistribution(normal_dist, normal_dist)
```

**Note**: 
!!! note
    Both `r_obs` and `err_resp` have the same type, eg. `MTResponse`. 


## Inference
We now have all the ingredients to perform inversion, except the sampler. This [page](https://turing.ml/dev/docs/using-turing/sampler-viz) provides a brief review of samplers `Turing.jl` provides. The number of posterior points to be sampled `n_samples`, the algorithm `mcmc_alg` and the distributions are brought together by `mcmc_cache`.
```julia
mcmc_alg = NUTS();
n_samples = 100;
mcache = mcmc_cache(modelD, respD, n_samples, mcmc_alg);
```

The posterior samples are then sampled by simply calling:
```julia
mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache, trans_utils = (m = log_tf,))
```
We know most of the variables here. `trans_utils` is the function that allows you to transform the samples drawn on `Uniform(-1, 5)` to on a log-scale on [$10^{-1}$, $10^5$]. Please check out the docstring of the function for more details.

## Copy-Pasteable code
```julia
using MT
using Distributions
using Turing
using LinearAlgebra


m_test = MTModel([100., 10., 1000.], [1e3, 1e3]);
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


modelD = MTModelDistribution(
    Product(
    [Uniform(-1., 5.) for i in eachindex(z)]
    ),
    vec(h)
);


n_samples = 50;
mcache = mcmc_cache(modelD, respD, 50, NUTS());

mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache, trans_utils = (m = log_tf,))
```

The obtained `mcmc_chain` contains the distributions that can be saved using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl).

```julia
using JLD2
JLD2.@save "file_path.jld2" mcmc_chains
```

**Note**: 
!!! note
    The returned chains will be sampled in the distribution specified by `modelD`. In the presented case, it will have values $\in [-1, 5]$ and we can get the values by `10. ^ value`.

The list of models can then be obtained from chains using
```
model_list = get_model_list(mcmc_chains, modelD, trans_utils = (m = log_tf,))
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

## Ideation
`max_depth` feature: This feature seems redundant because we're trying to force the interface between the last and the second last layer. If needed, one can simply specify so in the prior `modelDistribution` formulation. The last layer is a half-space and, therefore, doesn't require a thickness specified. 

`scale` feature: A better way to have this would be to specify in the distribution itself.

`responseDistribution`: Currently, it's not there but the way to go about this would be to have a similar `struct` as `modelDistribution`. We only need an example of a uniform distribution. Since sampling posterior using a `struct` is not possible, and our `forward` functions would return `response`s mostly, it's better to have the likelihood as a Normal Distribution. We'll open the option for covariance matrix rather than just the white noise diagonal matrix, and that itself will cover most of the requirements. 
We can drop a comment in the tutorial on how the `respDistribution` should be passed to be used as a likelihood function. 
In `Turing` model, we need to have the `responseDistribution` constructed around the `resp` evaluated on the sampled `model`. 


sampling `h`: The current implementation has `h` being sampled from the model distribution. More often than not, `h` is fixed throughout. While we can have tricks to specify a PDF to sample the same `h`, it'd be better to not sample it at all. We can have an option, for now but it's not going to be used a lot, especially when we look at 2D/ 3D, where grid points are fixed throughout.


For RTO, we can make a case and we have to actually look into how to best implement it but for general MCMC, we can have any misfit criteria and a `W` matrix. It's better to use `W` and if there is a loss function that doesn't require W, we'll simply not use it. In any case, it's an optional keyword that's just going to be used a lot more.

We'll have the `responseDistribution` and `modelDistribution` but these will be converted into `NamedTuple` before going through the mcmc sampler


Like we have `model_utils` to conviniently see if we want to infer on some of the model fields or all of them, a `resp_utils` won't make a lot of sense because 
* If you have the data, why not use it, and make it fix
* If you don't want to use that data, say `ϕ` for MT, simply exclude it in the `resp_fields`
* Now, whenever you want to have something like this, you always need to provide a `struct` containing all the fields for the response variable but you can fill the unwanted field with anything. Simple!
* What's the best way to include `resp_distribution`. It's fine to have it the way we have


Why we did not have Delta distribution?
* If something is constant, why sample it? Simply use the same value. Will be more efficient computationally, especially when the fixed variable is large, say `grid sizes` in 2D/3D MT forward computations.
* Also, it's more work


For most geophysical applications, we have our work done. But let's say when we want to extend the codes to eg., rock physics models, we'd simply need to extend the functions for making the `modelDistribution` and updating the `responseDistribution`. -->