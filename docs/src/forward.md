# Forward modeling

## Demo
Currently, only the recursion solution for MT forward modeling is supported. Once a model is defined, we can get the estimate as:

```@example forward_demo
using MT
using Plots
ρ= [500., 100., 400., 1000.];
h= [100., 1000., 10000.];
m= MTModel(ρ, h)

T= 10 .^(range(-1,5,length= 57));
ω= 2π./T;

resp= forward(m, ω);
nothing # hide
```

## Plots
Since MT `MTResponse` has two fields, and we also want to see the curves for other `MTModel`s, usually for inversion, it's easier to create a wrapper and then plot them.

```@example forward_demo
plt= prepare_plot(resp, ω, label = false)
plot_response(plt, margin= 5Plots.mm)
```

Another `MTResponse` can be overlain using:
```@example forward_demo
using MT
ρ= [100., 10., 400.];
h= [100., 10000.];
m2= MTModel(ρ, h)
resp2= forward(m2, ω)
prepare_plot!(plt, resp2, ω, label="2nd", plt_type = :plot)
plot_response(plt, margin= 5Plots.mm)
```

Note that we've plotted the second response curve as a line, which can simply be achieved by passing `plt_type = :plot`

## Benchmark

In-place operations for fast non-allocating computations are also supported. These would greatly speed up the inverse processes.

```@example forward_demo
using BenchmarkTools
@benchmark forward!(resp, m, ω)
```
```
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  25.875 μs … 119.209 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     28.209 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   28.593 μs ±   2.166 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

      ▂▄▅▅▅▆█▇▃▃▅▄▁▁▂▂▁                                        ▂
  ▂▂▄▆█████████████████▇▇██▇▇▇▇▇█▇█▆▅▆▆▆▅▅▅▄▄▂▃▄▄▄▄▅▅▄▄▅▃▄▄▅▄▅ █
  25.9 μs       Histogram: log(frequency) by time      38.9 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.
```

The above benchmark was done on Mac M1.