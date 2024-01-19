# # function forward!(ρₐ::T, m::T, ω::T, n::Union{Int64, Int32, Int16}) where T<:AbstractVector # m is the model vector
# function forward!(ρₐ, m, ω, n)
#     copyto!(ρₐ, (forward(model(m[1:n], m[n+1:end]), ω)).ρₐ);
#     # return resp.ρₐ
#     return nothing;
# end

# """
# Actually transpose(jacobian)
# J= Jacobian matrix, size: (length of forward vector, length of model vector)
# fn!= forward function
# m= model vector
# """
# function jacobian!(J::AbstractMatrix, m::AbstractVector, ω::AbstractVector, n::Union{Int64, Int32})
#     fn_tmp= zeros(size(J, 1)); # what about impedance?? let's worry about real number fns for now
#     # fn_m= zeros(size(J, 1));
#     for i in 1:length(m)
#         m[i]+= 1000000*eps(typeof(m[1]));
#         # @show typeof(m)<:AbstractVector
#         forward!(view(J, :, i), m, ω, n); # col-wise for performance
#         # @show view(J, :, i)
#         m[i]-= 2000000*eps(typeof(m[1]));
#         forward!(fn_tmp, m, ω, n);
#         broadcast!(-, view(J, :, i), view(J, :, i), fn_tmp);
#         # @show view(J, :, i)
        
#         rmul!(view(J, :, i), inv(2000000*eps(typeof(m[1]))));
#         m[i]+= 1000000*eps(typeof(m[1]));
#         # break;
#     end
#     return nothing;
# end

# # the following doesn't really help
# # linearized_forward(m, Js, m0, d0)= d0+ Js*(m.- m0);
# # function linearized_forward!(d, m, J, m0, d0)
# #     mul!(d, J, m- m0);
# #     broadcast!(+, d, d0, d);
# #     return nothing
# # end

# function gradient!(grad, m, dobs, Js, m0, d0, ω)
#     d= d0+ Js*(m.- m0); #linearized_forward(m, Js, m0, d0)
#     mul!(grad, transpose(Js), d- dobs)
#     rmul!(grad, inv(0.5* length(d)))
#     return nothing;
# end


# """
# linearize then least-square
# linearize the forward fn around m₀ (option for every iteration?, option for regularizer?)
# """
# function lls(dobs, m0, ω)
#     d0= forward(m0, ω).ρₐ
#     n= length(m0.ρ)
#     m0vec= [m0.ρ..., m0.h...];
#     Js= zeros(length(ω), length(m0vec));
#     jacobian!(Js, m0vec, ω, n);
    
#     J(mvec)= inv(length(d0))* norm(d0+ Js*(mvec.- m0vec).- dobs, 2);
#     grad!(grad, m)= gradient!(grad, m, dobs, Js, m0vec, d0, ω);
#     opt= optimize(J, grad!, m0vec, LBFGS(), inplace= true, Optim.Options(iterations= 20, show_trace= true))
#     return opt.minimizer
# end

# # m_init= model(m.ρ.*(1 .+ 0.1 .* randn(size(m.ρ))), m.h.*(1 .+ 0.1 .* randn(size(m.h))))

# # lls(resp.ρₐ, m_init, ω)