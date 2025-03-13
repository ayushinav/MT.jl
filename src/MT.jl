module MT
using LinearAlgebra
using Plots
using LinearSolve
using NonlinearSolve
using Optimization, OptimizationOptimJL
using Turing
using Distributions
using Statistics
using SpecialFunctions
using QuadGK
using UnPack
using ProgressMeter
import Base: show

include("abstract_types.jl")
include("probabilistic/init_distributions.jl")
include("models/mt.jl")
include("response/1dmt.jl")
include("forward/1dmt.jl")

include("rock_physics/conductivity/utils.jl")
include("rock_physics/conductivity/cache.jl")
include("rock_physics/conductivity/types.jl")
include("rock_physics/conductivity/forward.jl")

include("rock_physics/elastic/utils.jl")
include("rock_physics/elastic/cache.jl")
include("rock_physics/elastic/types.jl")
include("rock_physics/elastic/forward.jl")

include("rock_physics/viscous/utils.jl")
include("rock_physics/viscous/cache.jl")
include("rock_physics/viscous/types.jl")
include("rock_physics/viscous/forward.jl")

include("rock_physics/anelastic/utils.jl")
include("rock_physics/anelastic/cache.jl")
include("rock_physics/anelastic/types.jl")
include("rock_physics/anelastic/forward.jl")

include("rock_physics/mixing_phases.jl")
include("rock_physics/combine_models.jl")
include("utils.jl")
include("inverse/utils.jl")
include("inverse/bounds_transformation.jl")
include("inverse/jacobian.jl")
include("inverse/occam.jl")
include("inverse/inv.jl")
include("inverse/nl_inv.jl")
include("inverse/opt_inv.jl")
include("probabilistic/respDistribution.jl")
include("probabilistic/utils.jl")
include("probabilistic/inverse.jl")
include("probabilistic/rto.jl")
include("probabilistic/post_inv_utils.jl")
include("plots/utils.jl")
include("plots/plots.jl")
include("plots/prob_utils.jl")
include("rock_physics/pretty_printing.jl")

# export μ
export AbstractModel, AbstractResponse
export AbstractGeophyModel, AbstractGeophyResponse
export AbstractModelDistribution, AbstractResponseDistribution
export AbstractGeophyModelDistribution, AbstractGeophyResponseDistribution
export RockphyModelDistribution, RockphyResponseDistribution
export MTModel, MTResponse
export RockphyCond
export get_Z, get_appres, get_phase, forward!, forward
# export zero, copy
export plot_response, prepare_plot, prepare_plot!
export plot_model, plot_model!
export sigmoid, d_sigmoid, inverse_sigmoid, transform_utils, default_tf, log_tf, lin_tf
export mt_jacobian_cache, jacobian_mt, jacobian!
export occam_cache, Occam, nl_cache, NonlinearAlg, opt_cache, OptAlg, linsolve!, occam_step!
export inverse!
export ∂, χ², linear_utils, inverse_utils
export normal_dist, uniform_dist
export MTModelDistribution, MTResponseDistribution
export RockphyModelDistribution, RockphyResponseDistribution
export SEO3, UHO2014, Jones2012, Poe2010, Yoshino2009, Wang2006, const_matrix
export Ni2011, Sifre2014, Gaillard2008
export anharmonic, anharmonic_poro, SLB2005
export HZK2011, HK2003, xfit_premelt
export eburgers_psp, andrade_psp, andrade_analytical, premelt_anelastic, xfit_mxw
export model_multiphase, consturct_model_multiphase
export construct_mixing_models, mixing_models, HS1962_plus, HS1962_minus, single_phase, MAL
export construct_model_multi_rp, model_multi_rp
export mcmc_cache, rto_cache
export stochastic_inverse, get_model_list
export pre_image, get_kde_image, get_mean_std_image

end
