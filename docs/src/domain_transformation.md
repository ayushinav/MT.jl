# Domain transformation

## Introduction

Throughout the framework, you will come across `transform_utils`, whether be it deterministic inverse, or probabilistic. Even in simple forward calculations, we have them embedded.

Firstly, let's understand its structure with the example of `log_tf`. A `transform_utils` type variable, say `trans_util` has 4 parameters, the second of which, referred by `tf` is used to move from one space to another. For consistency, we call this forward transformation, that is, it will move the variables to the log-space (on the base 10).

```julia
x = 10.0
log_tf.tf(x) ## = 1.
```

Then we have inverse transformation given by `itf`, and it will allow you to get the values back.

```julia
x = 2.0
log_tf.tf(x) ## = 100.
```

We, then also have `dtf`, which computes the derivative of the forward tramnsformation

```julia
x = 10.0
log_tf.tf(x) ## 
```

We also have `p` which is a vector to parameterize these functions if needed, eg., `sigmoid_tf` has `p` which define the lower and upper bounds of the scaled sigmoid. This is used in bounding variables without putting an explicit constraint on them.

## Model transformation

Model transformation is where these transformations really help. As an eg., for MT inversion, we want our bounds of log resistivity to be between `-1` and `3`. We then define a new `transform_utils` variable as:

```julia
my_trans_util = transform_utils([-1.0, 3.0], sigmoid, inverse_sigmoid, d_sigmoid);
```

and then pass this in the `model_transform_utils` parameter of the functions such as [`inverse!`](@ref), [`stochastic_inverse`](@ref) and many others.

Explain that this should mostly be used for bounding, example of regularization.

## Response transformation

Similar to `model_trans_utils`, a lot of functions, including `forward` have the option to pass in `response_trans_utils`. The working is fairly similar in the aspect that we transform the data. This can be useful in inversion schemes when we have data in different varying in different domains, eg. apparent resistivity varying smoothly on the log-scale and phase on the linear scale.

Explain how data to be given should already be in the transformed space
