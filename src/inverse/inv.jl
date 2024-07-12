"""
    function inverse!(mₖ::model,
            robs::response,
            vars::Vector{Float64},
            alg_cache::occam_cache;
            W= nothing, # Weight matrix
            max_iters= 20, χ2=1.,
            response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(robs))],
            model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(m₀))], # this will not be used but for the sake of generality for all inverse algs
            trans_utils::transform_utils= sigmoid_tf,
            verbose= true
        ):
updates `mₖ` using occam iteration to fit `robs` within a misfit of `χ2`, by default set to 1.0. 
### Variables:
* `mₖ::model`: Inital model guess, will be updated during the inverse process
* `robs::response`: response to invert for
* `vars::Vector{Float64}`: variables required for forward modeling, eg., `ω` for MT
* `alg_cache::occam_cache`: deterimines the algorithm to be performed for inversion
* `W= nothing`: Weight matrix, will be `I` if nothing is provided
* `max_iters= 20`: maximum number of iterations
* `χ2=1.`: threshold misfit
* `response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(robs))]`: choose data of response to perform inversion on, eg., ρₐ for MT, by default chooses all the data (ρₐ and ϕ)
* `model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(m₀))]`: will generally be fixed, see docs for details
* `trans_utils::transform_utils= sigmoid_tf`: for bounds transformation,
* `verbose`: whether to print updates after each iteration, defaults to true
* `mᵣ` : model to be regularized against

### Returns:
return message in the form of `return_code` and updates `mₖ` in-place.

### Example:
`inverse!(m_occam, r_obs, Occam([1e-2, 1e6]))`
"""
function inverse!(mₖ::model1,
        robs::response,
        vars::Vector{Float64},
        alg_cache::occam_cache;
        W= nothing, # Weight matrix
        max_iters= 20, χ2=1.,
        response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(robs))],
        model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(mₖ))], # this will not be used but for the sake of generality for all inverse algs
        trans_utils::transform_utils= sigmoid_tf,
        verbose::Bool= true,
        mᵣ::model2 = zero(mₖ)
    ) where {model1 <: AbstractGeophyModel, model2 <: AbstractGeophyModel, response <: AbstractGeophyResponse}

    prec= eltype(mₖ.m);
    model_fields= [:m];


    n_model= length(mₖ.m); 
    n_vars= length(vars);
    n_resp= length(response_fields)* n_vars;

    (W == nothing) && (W= prec.(I(n_resp)));

    lin_utils= linear_utils(view(mₖ.m, :), zeros(prec, n_resp), zeros(prec, n_resp, n_model));

    respₖ= zero_abstract(robs); # MTResponse{AbstractVector{prec}}([zero(vars) for k in fieldnames(typeof(robs))]...)
    jₖ= jacobian_mt(fieldnames(typeof(robs)), eltype(vars));

    for (i,k) ∈ enumerate(response_fields)
        setfield!(jₖ, k, view(lin_utils.Jₖ, (i-1)*n_vars+ 1:i*n_vars, :));
        setfield!(respₖ, k, view(lin_utils.Fₖ, (i-1)*n_vars+1: i*n_vars));
    end
    
    mtjc= mt_jacobian_cache(vars);

    inv_utils= inverse_utils(∂(n_model), W, reduce(vcat, [copy(getfield(robs, k)) for k ∈ response_fields]));

    mₖ₊₁= copy(mₖ); # MTModel([copy(getfield(mₖ, k)) for k ∈ fieldnames(typeof(mₖ))]...)
    respₖ₊₁= copy(respₖ); # MTResponse([copy(getfield(respₖ, k)) for k ∈ fieldnames(typeof(respₖ))]...)

    lin_prob= LinearProblem(inv_utils.D'inv_utils.D,
        lin_utils.Jₖ'*(inv_utils.dobs+ lin_utils.Jₖ*lin_utils.mₖ- lin_utils.Fₖ));
    linsolve_prob= init(lin_prob, assumptions=LinearSolve.OperatorAssumptions(true, condition= LinearSolve.OperatorCondition.WellConditioned));

    forward!(respₖ, mₖ, vars); # for the first iteration
    itr= 1;
    chi2= prec(1e6);

    μ_last = 0.;
    while itr<= max_iters
        verbose && (print("$itr: ");)
        jacobian!(jₖ, mₖ, vars, mtjc, model_fields= model_fields, response_fields= response_fields);
        copyto!(lin_utils.Jₖ, lin_utils.Jₖ.*lin_utils.mₖ'.* log(10));
        for k in model_fields # to computational domain
            getfield(mₖ, k).= trans_utils.itf.(log10.(getfield(mₖ, k)));
        end

        μ_last = occam_step!(mₖ₊₁, # to store the next update, which will eventually be copied to mₖ
            mᵣ, # model to be regularized against
            respₖ₊₁, # to store the response for mₖ₊₁, for error calculation and anything
            vars, # to compute the forward model
            χ2, # threshold chi-squared error that needs to be met
            alg_cache.μgrid, # for gridsearch of μ for Occam
            lin_utils, # contains the mₖ, Jₖ, Fₖ associate with the current iteration
            inv_utils, # contains D= ∂(n), W and dobs
            trans_utils, # to  transform to and from the computational domain
            linsolve_prob, # for faster inverse operations
            model_fields= model_fields,
            response_fields= response_fields,
            verbose= verbose, 
        )

        for k in model_fields # copying things to mₖ
            getfield(mₖ, k).= getfield(mₖ₊₁, k)
        end
        forward!(respₖ, mₖ, vars);
        chi2= χ²(reduce(vcat, [copy(getfield(respₖ, k)) for k ∈ response_fields]), inv_utils.dobs, W= inv_utils.W);
        if chi2< χ2 break; end
        itr+=1;
    end

    return return_code(
        chi2<=χ2,
        (μ = μ_last,),
        mₖ,
        χ2,
        chi2
    )    
end