# Probabilistic inversion

## Ideation
`max_depth` feature: This feature seems redundant because we're trying to force the interface between the last and the second last layer. If needed, one can simply specify so in the prior `modelDistribution` formulation. The last layer is a half-space and, therefore, doesn't require a thickness specified. 

`scale` feature: A better way to have this would be to specify in the distribution itself.

`responseDistribution`: Currently, it's not there but the way to go about this would be to have a similar `struct` as `modelDistribution`. We only need an example of a uniform distribution. Since sampling posterior using a `struct` is not possible, and our `forward` functions would return `response`s mostly, it's better to have the likelihood as a Normal Distribution. We'll open the option for covariance matrix rather than just the white noise diagonal matrix, and that itself will cover most of the requirements. 
We can drop a comment in the tutorial on how the `respDistribution` should be passed to be used as a likelihood function.


sampling `h`: The current implementation has `h` being sampled from the model distribution. More often than not, `h` is fixed throughout. While we can have tricks to specify a PDF to sample the same `h`, it'd be better to not sample it at all.