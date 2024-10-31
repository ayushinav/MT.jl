
# TODO : change these to f-5 and f-19
const global boltz_k = 8.617e-5;
const global charge_e = 1.602e-19;

"""
T is in K
"""
function forward(m::model, p) where {model <: SEO3}
    fO₂ = 10. ^(-24441.9/m.T + 13.296)

    # sT 
    bfe = 5.06e24 * exp(-0.357* inv(boltz_k*m.T));
    bmg = 4.58e26 * exp(-0.752* inv(boltz_k*m.T));
    ufe = 12.2e-6 * exp(-1.05* inv(boltz_k*m.T));
    umg = 2.72e-6 * exp(-1.09* inv(boltz_k*m.T));
    concFe = bfe + 3.33e24 * exp(-0.02* inv(boltz_k*m.T)) * fO₂ ^(1/6); 
    concMg = bmg + 6.21e30*exp(-1.83* inv(boltz_k*m.T)) * fO₂ ^(1/6); 
    σ = concFe *ufe * charge_e + 2. * concMg * umg * charge_e;

    return σ

end

"""
T is in K
"""
function forward(m::Ni2011, p)

    Tcorr = 1146.8 # K
    D = 0.006 # unitless, Partition coefficient {ol/melt}

    ls = 2.172 - (860.82 - 204.46 * sqrt(m.Ch2o_m/1e4)) * inv(m.T - Tcorr);
    σ = 10. ^ ls

    return σ

end
