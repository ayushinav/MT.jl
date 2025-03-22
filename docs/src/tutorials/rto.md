## RTO-TKO

RTO-TKO is a stochastic framework introduced in the MT community by Blatter et al., 2022 ((a)[ https://doi.org/10.1093/gji/ggac241] and (b)[https://doi.org/10.1093/gji/ggac242])

Similar to fixed discretization scheme, the grid sizes do not change. RTO-TKO chooses a perturbation in the prior space and optimzes for a response sampled from the response domain.

Let's have

```math
m = [m_1, m_2, m_3, ... , m_N] \\
h = [h_1, h_2, h_3, ... , h_{N_1}]
```

and

```math
d = \mathcal{F}(m)
```

Let ``C_d`` be the data covariance matrix and we want to explore uncertainty for the following misfit function:

```math
J(m) = [\mathcal{F}(m) - d]^T C_d [\mathcal{F}(m) - d] + \mu (Lm)^T(Lm)
```

where ``\mu`` is the regularization weight and ``L`` is the derivative matrix. RTO-TKO explores the uncertainty in ``\mu``-space instead of fixing it and gives a family of models that fit the data. The algorithm works as:

Solving for

```math
\begin{align*}
& J(m) = [\mathcal{F}(m) - d]^T C_d [\mathcal{F}(m) - d] + \mu (Lm)^T(Lm) \\

& 1) \quad \text{Solve for } \; m_{i+1} \\
& \text{Sample  } \; \tilde{d} \sim \mathcal{N}(d, C_d) \text{ and } \tilde{m} \sim \mathcal{N}(0, \frac{1}{\mu}(L^T L))\\

& \text{Solve} \\
& m_{i+1} = \argmin_{m_{i+1}} \quad [\mathcal{F}(m_{i+1}) - \tilde{d}]^T C_d [\mathcal{F}(m_{i+1}) - \tilde{d}] + \mu_i [L(m_{i+1} - \tilde{m})]^T[L(m_{i+1} - \tilde{m})] \\

& 2) \quad \text{Solve for } \; \mu_{i+1} \\
& \text{Sample  } \; \tilde{d} \sim \mathcal{N}(d, C_d) \\

& \text{Solve} \\
& \mu_{i+1} = \argmin_{\mu_{i+1}} \quad [\mathcal{F}(m_{i+1}) - \tilde{d}]^T C_d [\mathcal{F}(m_{i+1}) - \tilde{d}] - log(p(\mu_{i+1})) \\
\end{align*}
```

!!! note
    
      - Usually, the prior of ``\mu``is a uniform distribution and we do not have to compute the log pdf term
      - Implementing the above from scratch might not be trivial because of ``L^T L`` being non-invertible

## Copy-Pasteable code

```julia
using MT
using Distributions
using Turing
using LinearAlgebra

m_test = MTModel(log10.([100.0, 10.0, 1000.0]), [1e3, 1e3]);
f = 10 .^ range(-4; stop=1, length=25);
ω = vec(2π .* f);

r_obs = forward(m_test, ω);

err_phi = asin(0.01) * 180 / π .* ones(length(ω));
err_appres = 0.02 * r_obs.ρₐ;
err_resp = MTResponse(err_appres, err_phi);

r_obs.ρₐ .= r_obs.ρₐ .+ err_appres;
r_obs.ϕ .= r_obs.ϕ .+ err_phi;

respD = MTResponseDistribution(normal_dist, normal_dist);

z = 10 .^ collect(range(1; stop=4, length=100));
h = diff(z);

# fixed discretization
modelD = MTModelDistribution(
    Product(
        [Uniform(-1.0, 5.0) for i in eachindex(z)]
    ),
    vec(h)
);

n_samples = 50;
r_cache = MT.rto_cache(m_rto, [1e-6, 1e2], Occam(), 50, 1000, 1.0, [:ρₐ, :ϕ], false)

rto_chain = stochastic_inverse(r_obs, err_resp, ω, r_cache)

mt_chain = Turing.Chains(
    (rto_chains.value.data[:, 1:length(z), :]), [Symbol("ρ[$i]") for i in 1:length(z)]);
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
model_list = get_model_list(mt_chains, modelD)
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
kde_img = get_kde_image(pre_img..., false; xscale=:identity, yscale=:identity, yflip=true)
```

We can also obtain the mean and 1 std deviation bounds as:

```julia
mean_std_plt_lin = get_mean_std_image(pre_img...; yscale=:identity)
```
