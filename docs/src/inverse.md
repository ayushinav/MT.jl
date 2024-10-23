<!-- # Inversion

## Brief introduction
Inverse problems in geophysics are notoriously ill-posed with non-unique solutions. MT inversion is no different. Occam inversion is one of the most popular inversion algorithms to approach MT problems, and it comes with its own challenges. Not much to argue when dealing with compiling fortran codes and even less when one realizes that the Occam code doesn't work very well on Mac chips (not sure about AMD ones). Here we have the Occam inversion code on a test synthetic dataset, as easy as it can be. 

## Demo
We start with defining models:

```@example inverse_demo
using MT

ρ= [100., 10., 400., 1000.]
h= [100., 1000., 10000.]
m= MTModel(ρ, h)

T= 10 .^(range(-1,5,length= 57))
ω= 2π./T

plot_model(m)
```

and getting data.

```@example inverse_demo
resp= forward(m, ω)

plt= prepare_plot(resp, ω, label= false)
plot_response(plt)
```

One can add errors and everything but here we just put up a simple example.
```@example inverse_demo
h_test= 10 .^range(0., 5., length= 50)
ρ_test= 5e2 .*ones(length(h_test)+1)

m_test= MTModel(ρ_test, h_test);
nothing # hide
```

Now showtime, summon Occam.

```@example inverse_demo
inverse!(m_test, resp, ω, MT.Occam(μgrid = [1e-2, 1e6]), max_iters= 50)
nothing # hide
```

Let's look at how our model looks like:
```@example inverse_demo
plt= plot_model(m_test, label= "inverted");
plot_model!(plt, m, label= "true")
```

and how well it fits the data:
```@example inverse_demo
resp_true= forward(m, ω)
resp_test= forward(m_test, ω);

plt= prepare_plot(resp_true, ω, label= "true", plt_type = :scatter)
prepare_plot!(plt, resp_test, ω, label= "inverted", plt_type = :plot, linewidth= 3)
plot_response(plt,legend=:bottomright)
```

## Benchmark

And let's finish up with a benchmark on Mac M1. We run 10 iterations to just get an estimate on how things run. Benchmarking is not possible using `@benchmark` because the consecutive operations will converge in a single iteration because of the in-place updates.

```@example inverse_demo
h_test= 10 .^range(0., 5., length= 50);
ρ_test= 5e2 .*ones(length(h_test)+1);

m_test= MTModel(ρ_test, h_test);
@time inverse!(m_test, resp, ω, Occam(), max_iters= 10)
nothing # hide
```

Finally, the `inverse!` returns a `retcode` that might be beneficial while building customized codes for those specific cases.
```@example inverse_demo
retcode = inverse!(m_test, resp, ω, Occam(), max_iters= 10)
nothing # hide
```

 That easy and that fast! -->