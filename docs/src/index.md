# `MT.jl`

`MT.jl` is supposed to be a high performance code for doing forward and inverse modeling in julia. We hope to write the code structure such that any other geophysical survey can also be used and we can tend towards a joint forward and inverse modeling library.

That being said, we currently support 1D forward modeling and Occam inversion only.


## Forward modeling

While forward modeling typically requires solving a PDE obtained using the quasi-static approximation, in 1D, we are fortunate to have the solution for surface impedance in a more analytical form. Currently, this is what is supported.

### Ideation
Until this will go public, we can also use this page to develop the ideas for code structure for forward modeling. We know that in future (not so far future), we will have 2D (atleast) forward model. It might probably make sense to have a different folder for `forward` codes, and then have FD, FE or analytical soln in there for 1D, and/ or 2D.


## Inverse modeling
No surprises here that we are almost always trying to solve for an under-determined system. The most popular framework in 1D is Occam and this is what is currently supported at the moment.

### Ideation
 We can simply add gradient descent (though not useful) but hopefully, we can add more neural network oriented framework for our approach. Impyling, occam codes go in a sub folder in the `inverse` codes. This can be better approached using `NonLinearSolve.jl`.

Check out the implementation details:
* [model](model.md)
* [forward](forward.md)
* [inverse](inverse.md)
<!-- * [probabilistic_inverse](probabilisitc_inverse.md) -->
* [interface_guide](interface_guide.md)
* [rock physics](rock_physics.md)

<!-- ```@meta
CurrentModule = MT
```

# MT

Documentation for [MT](https://github.com/ayushinav/MT.jl).

```@index
```

```@autodocs
Modules = [MT]
``` -->
