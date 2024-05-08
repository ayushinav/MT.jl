# Probabilistic inversion

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


For most geophysical applications, we have our work done. But let's say when we want to extend the codes to eg., rock physics models, we'd simply need to extend the functions for making the `modelDistribution` and updating the `responseDistribution`.