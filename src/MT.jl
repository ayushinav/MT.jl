module MT
using LinearAlgebra
using Plots
using LinearSolve
using Turing
using Distributions

include("models/1dmt.jl")
include("response/1dmt.jl")
include("forward/1dmt.jl")
include("plots/plots.jl")
include("inverse/utils.jl");
include("inverse/bounds_transformation.jl");
include("inverse/jacobian.jl");
include("inverse/occam.jl");
include("inverse/inv.jl");
include("probabilistic/init_distributions.jl")
include("probabilistic/gaussian_distribution.jl")
include("probabilistic/utils.jl")
include("probabilistic/inverse.jl")


# export μ
export model, response
export get_Z, get_appres, get_phase, forward!, forward
export plot_response, prepare_plot, prepare_plot!
export plot_model, plot_model!
export sigmoid, d_sigmoid, inverse_sigmoid, transform_utils, default_tf, log_tf
export mt_jacobian_cache, jacobian_mt, jacobian!
export occam_cache, Occam, linsolve!, occam_step!
export ∂, χ², linear_utils, inverse_utils;
export gaussian_responseDistribution, update_responseDistribution!
export modelDistribution, responseDistribution
export mcmc_cache
export inverse!; # both for occam and mcmc

end
