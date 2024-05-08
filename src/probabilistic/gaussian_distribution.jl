function gaussian_responseDistribution(resp::response, err_resp::response{T}; resp_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(resp))]) where T <: AbstractVector

    arr = [];
    for k ∈ fieldnames(typeof(resp))
        if k in resp_fields
            push!(arr, MultivariateNormal(
            copy(getfield(resp, k)), copy(Diagonal(getfield(err_resp, k)))
            ))
        else
            push!(arr, nothing)
        end
    end

    return responseDistribution(arr...)
end


function gaussian_responseDistribution(resp::response, err_resp::response{T}) where T <: AbstractMatrix
    responseDistribution(
        ([MultivariateNormal(
        copy(getfield(resp, k)), copy(getfield(err_resp, k))
        ) 
        k ∈ fieldnames(typeof(resp))]...)
    )
end

function update_responseDistribution!(respD::responseDistribution{T}, resp::response;
    response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(resp))]) where T <: Union{MultivariateNormal, Nothing}

    for k in response_fields
        getfield(getfield(respD, k), :μ) .= copy(getfield(resp, k))
    end
end