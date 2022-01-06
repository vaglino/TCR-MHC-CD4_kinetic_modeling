# analytical solution

function two_state_catch_bond_analytical_w_cons(p,M,t;f=0.0)

    x₁₀⁰,      # k₋ᵒf,x₋ᵒf
    k₂₀⁰,x₂₀⁰,      # k₋ᵒs,x₋ᵒs
    k₁₂⁰,x₁₂⁰,      # kᵒa,xᵒa,
    k₂₁⁰,x₂₁⁰ = p   # k₋ᵒa,x₋ᵒa
    # @show p
    # @show M.cons
    k₁₀⁰ = M.cons[1]
    # @show cons[1]

    k₁₀ = bell_diss_kT(k₁₀⁰, x₁₀⁰, f) # fast diss form L-R
    k₂₀ = bell_diss_kT(k₂₀⁰, x₂₀⁰, f) # slow diss from L-R*
    k₁₂ = bell_diss_kT(k₁₂⁰, x₁₂⁰, f)   # activation from L-R to L-R*
    k₂₁ = bell_ass_kT( k₂₁⁰, x₂₁⁰, f)  # deactivation from L-R* to L-R

    # initial states
    # B₁⁰ =(k₂₁⁰ * k₁₀⁰) / (k₂₁⁰ * k₁₀⁰ + k₁₂⁰ * k₂₀⁰)
    # B₂⁰ = 1 - B₁⁰
    # println("old")
    # @show B₁⁰_,B₂⁰_
    B₁⁰,B₂⁰ = initial_u(p,M)
    # if M.u₀_type == "equilibrium"
    #     B₁⁰,B₂⁰ = initial_state(k₁₂⁰,k₂₁⁰)
    # elseif M.u₀_type == "detailed_balance"
    #     B₁⁰,B₂⁰ = initial_state(k₁₀⁰,k₂₀⁰,k₁₂⁰,k₂₁⁰)
    #     # println("new")
    #     # @show B₁⁰,B₂⁰
    # end
    
    b = k₂₁ + k₂₀ + k₁₂ + k₁₀
    c = k₂₁*k₁₀ + k₁₀*k₂₀ + k₁₂*k₂₀
    λ₁ = λ1(b,c)
    λ₂ = λ2(b,c)
    C₁ = (k₂₁ + k₁₂ + B₁⁰*k₂₀ + B₂⁰*k₁₀ - λ₁) / (λ₂ - λ₁)
    C₂ = 1. - C₁

    pdf = C₁*λ₁*exp(-λ₁*t) + C₂*λ₂*exp(-λ₂*t)
    # pdf >= 0. || @show B₁⁰,B₂⁰
    return pdf
end

λ1(b,c) = (b + sqrt(b^2 - 4*c)) / 2
λ2(b,c) = (b - sqrt(b^2 - 4*c)) / 2


function two_state_catch_bond_analytical(p,M,t;f=0.0)

    k₁₀⁰,x₁₀⁰,      # k₋ᵒf,x₋ᵒf
    k₂₀⁰,x₂₀⁰,      # k₋ᵒs,x₋ᵒs
    k₁₂⁰,x₁₂⁰,      # kᵒa,xᵒa,
    k₂₁⁰,x₂₁⁰ = p   # k₋ᵒa,x₋ᵒa

    # @show cons[1]

    k₁₀ = bell_diss_kT(k₁₀⁰, x₁₀⁰, f) # fast diss form L-R
    k₂₀ = bell_diss_kT(k₂₀⁰, x₂₀⁰, f) # slow diss from L-R*
    k₁₂ = bell_diss_kT(k₁₂⁰, x₁₂⁰, f)   # activation from L-R to L-R*
    k₂₁ = bell_ass_kT( k₂₁⁰, x₂₁⁰, f)  # deactivation from L-R* to L-R

    B₁⁰,B₂⁰ = initial_u(p,M)
    # @show B₁⁰,B₂⁰

    b = k₂₁ + k₂₀ + k₁₂ + k₁₀
    c = k₂₁*k₁₀ + k₁₀*k₂₀ + k₁₂*k₂₀
    λ₁ = λ1(b,c)
    λ₂ = λ2(b,c)
    C₁ = (k₂₁ + k₁₂ + B₁⁰*k₂₀ + B₂⁰*k₁₀ - λ₁) / (λ₂ - λ₁)
    C₂ = 1 - C₁

    pdf = C₁*λ₁*exp(-λ₁*t) + C₂*λ₂*exp(-λ₂*t)
end

function two_states_analytical(p,cons,t;f=0.0)
    k₁₀⁰,x₁₀⁰,      # k₋ᵒf,x₋ᵒf
    k₂₀⁰,x₂₀⁰ = p      # k₋ᵒs,x₋ᵒs

    k₁₀ = bell_diss_kT(k₁₀⁰, x₁₀⁰, f) # fast diss form L-R
    k₂₀ = bell_diss_kT(k₂₀⁰, x₂₀⁰, f) # slow diss from L-R*

    pdf = 0.5 * k₁₀*exp(-k₁₀*t) + 0.5 * k₂₀*exp(-k₂₀*t)
end

function trimolecular_indep_catch_analytical_w_cons(p,M,t;f=0.0)
    # @show M.cons
    k_bi₁ = M.cons[1:8]
    k_bi₂ = M.cons[9:16]
    k₋ᵒf₃ = M.cons[17]
    M_bi₁ = deepcopy(M);  M_bi₁.cons = []; M_bi₁.u₀ = M.u₀[1:2]; M_bi₁.u₀_type = []
    M_bi₂ = deepcopy(M);  M_bi₂.cons = []; M_bi₂.u₀ = M.u₀[3:4]; M_bi₂.u₀_type = []
    M₃ = deepcopy(M);     M₃.cons = k₋ᵒf₃;
    if M.u₀_type == []
        M₃.u₀ = M.u₀[5:6]
    else
        M₃.u₀ = []; # and u₀_type says same as specified in M
    end

    pdf₁ = two_state_catch_bond_analytical(k_bi₁,M_bi₁,t;f=f)
    pdf₂ = two_state_catch_bond_analytical(k_bi₂,M_bi₂,t;f=f)
    pdf₃ = two_state_catch_bond_analytical_w_cons(p,M₃,t;f=f)

    # pdf = pdf₁ + pdf₂ + pdf₃
    
    pdf = u₀_mix[1]* pdf₁ + u₀_mix[2]*pdf₂ + u₀_mix[3]*pdf₃
    # pdf = 0.454* pdf₁ + 0.290*pdf₂ + 0.256*pdf₃
    return pdf
end

function trimolecular_indep_catch_analytical(p,M,t;f=0.0)
    # @show M.cons
    k_bi₁ = M.cons[1:8]
    k_bi₂ = M.cons[9:16]
    # k₋ᵒf₃ = M.cons[17]
    M_bi₁ = deepcopy(M);  M_bi₁.cons = []; M_bi₁.u₀ = M.u₀[1:2]; M_bi₁.u₀_type = []
    M_bi₂ = deepcopy(M);  M_bi₂.cons = []; M_bi₂.u₀ = M.u₀[3:4]; M_bi₂.u₀_type = []
    M₃ = deepcopy(M);
    if M.u₀_type == []
        M₃.u₀ = M.u₀[5:6]
    else
        M₃.u₀ = []; # and u₀_type says same as specified in M
    end

    pdf₁ = two_state_catch_bond_analytical(k_bi₁,M_bi₁,t;f=f)
    pdf₂ = two_state_catch_bond_analytical(k_bi₂,M_bi₂,t;f=f)
    pdf₃ = two_state_catch_bond_analytical(p,M₃,t;f=f)

    # pdf = pdf₁ + pdf₂ + pdf₃
    
    pdf = u₀_mix[1]* pdf₁ + u₀_mix[2]*pdf₂ + u₀_mix[3]*pdf₃
    # pdf = 0.454* pdf₁ + 0.290*pdf₂ + 0.256*pdf₃
    return pdf
end


function slip_bond_analytical(p,M,t;f=0.0)
    k₋₁⁰,x₋₁⁰ = p   
    k₋₁ = bell_diss_kT(k₋₁⁰,x₋₁⁰, f) # diss from L-R
  
    pdf = k₋₁*exp(-k₋₁*t)
end

function trimolecular_2_species_indep_catch_analytical(p,M,t;f=0.0)
    # @show M.cons
    k_bi = M.cons[1:8]


    M_bi = deepcopy(M);  M_bi.cons = []; M_bi.u₀ = M.u₀[1:2]; M_bi.u₀_type = []
    M₃ = deepcopy(M);
    if M.u₀_type == []
        M₃.u₀ = M.u₀[3:4]
    else
        M₃.u₀ = []; # and u₀_type says same as specified in M
    end

    pdf_bi = two_state_catch_bond_analytical(k_bi,M_bi,t;f=f)
    pdf₃ = two_state_catch_bond_analytical(p,M₃,t;f=f)

    # pdf = pdf₁ + pdf₂ + pdf₃
    
    pdf = u₀_mix[1]*pdf_bi + u₀_mix[2]*pdf₃
    return pdf
end

# function trimolecular_catch_analytical(p,M,t;f=0.0)
#     k₋ᵒf,x₋ᵒf,
#         k₋ᵒs,x₋ᵒs,
#         kᵒa,xᵒa,
#         k₋ᵒa,x₋ᵒa = cons
#     k₂⁰,x₂⁰,  
#         k₋₂⁰,x₋₂⁰ = p 

#     kf = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
#     ks = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
#     ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
#     k_a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)  # deactivation from L-R* to L-R
#     k2  = bell_ass_kT(k₂⁰, x₂⁰, f)   # X association w/ L-R
#     k_2 = bell_diss_kT(k₋₂⁰, x₋₂⁰, f) # X dissociation from L-R-X
#     # model

#     R = (exp(-t*(k_a + ks))*(k2*k_a*uA - ks^2*uA - 2*k_a^2*uA + k2*ks*uA + 2*k_2*k_a*uA - k_2*ka*uA + k_a*ka*uA - k_2*kf*uA + k_a*kf*uA - k_2*ka*uR - k_2*ka*uRX + k_a*ka*uR + k_2*ks*uA - 3*k_a*ks*uA + ka*ks*uA + kf*ks*uA + ka*ks*uR))/(k2*k_a + k2*ks + 2*k_2*k_a - k_2*ka + k_a*ka - k_2*kf + k_a*kf + k_2*ks - 3*k_a*ks + ka*ks + kf*ks - 2*k_a^2 - ks^2) + (exp(-(t*(k2 + k_2 - k_a + ka + kf - (k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/2)*((k_2*ka*(k2 + k_2 - k_a - ks))/(k2*(k2*k_a + k2*ks + 2*k_2*k_a - k_2*ka + k_a*ka - k_2*kf + k_a*kf + k_2*ks - 3*k_a*ks + ka*ks + kf*ks - 2*k_a^2 - ks^2)) + (ka*(k_a - k_2 + ks)*(k2/2 + k_2/2 - k_a/2 + ka/2 + kf/2 - (k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)/2))/(k2*(k2*k_a + k2*ks + 2*k_2*k_a - k_2*ka + k_a*ka - k_2*kf + k_a*kf + k_2*ks - 3*k_a*ks + ka*ks + kf*ks - 2*k_a^2 - ks^2)))*(2*k2*uR + k2*uRX - k_2*uRX - k_a*uRX + ka*uRX + kf*uRX + uRX*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/(2*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)) - (exp(-(t*(k2 + k_2 - k_a + ka + kf + (k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/2)*((k_2*ka*(k2 + k_2 - k_a - ks))/(k2*(k2*k_a + k2*ks + 2*k_2*k_a - k_2*ka + k_a*ka - k_2*kf + k_a*kf + k_2*ks - 3*k_a*ks + ka*ks + kf*ks - 2*k_a^2 - ks^2)) + (ka*(k_a - k_2 + ks)*(k2/2 + k_2/2 - k_a/2 + ka/2 + kf/2 + (k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)/2))/(k2*(k2*k_a + k2*ks + 2*k_2*k_a - k_2*ka + k_a*ka - k_2*kf + k_a*kf + k_2*ks - 3*k_a*ks + ka*ks + kf*ks - 2*k_a^2 - ks^2)))*(2*k2*uR + k2*uRX - k_2*uRX - k_a*uRX + ka*uRX + kf*uRX - uRX*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/(2*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2))
    
#     A = (exp(-(t*(k2 + k_2 - k_a + ka + kf + (k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/2)*(2*k2^2*uR + 2*k2*uR*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2) - 2*k2*k_2*uR - 4*k2*k_2*uRX - 2*k2*k_a*uR + 2*k2*ka*uR + 2*k2*kf*uR))/(4*k2*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)) + (exp(-(t*(k2 + k_2 - k_a + ka + kf - (k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/2)*(2*k2*uR*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2) - 2*k2^2*uR + 2*k2*k_2*uR + 4*k2*k_2*uRX + 2*k2*k_a*uR - 2*k2*ka*uR - 2*k2*kf*uR))/(4*k2*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2))
    
#     RX = (exp(-(t*(k2 + k_2 - k_a + ka + kf - (k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/2)*(2*k2*uR + k2*uRX - k_2*uRX - k_a*uRX + ka*uRX + kf*uRX + uRX*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/(2*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)) - (exp(-(t*(k2 + k_2 - k_a + ka + kf + (k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/2)*(2*k2*uR + k2*uRX - k_2*uRX - k_a*uRX + ka*uRX + kf*uRX - uRX*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2)))/(2*(k2^2 + 2*k2*k_2 - 2*k2*k_a + 2*k2*ka + 2*k2*kf + k_2^2 + 2*k_2*k_a - 2*k_2*ka - 2*k_2*kf + k_a^2 - 2*k_a*ka - 2*k_a*kf + ka^2 + 2*ka*kf + kf^2)^(1/2))
    


# end


# uR*exp(-t*(ka - k_a + kf))*(ka - k_a + kf) - (exp(-t*(k_a + ks))*(k_a + ks)*(ka*uA - 2*k_a*uA + kf*uA + ka*uR - ks*uA))/(2*k_a - ka - kf + ks) + (ka*uR*exp(-t*(ka - k_a + kf))*(ka - k_a + kf))/(2*k_a - ka - kf + ks)

function two_state_catch_bond_ana_matlab(p,M,t;f=0.0)

    k₁₀⁰,x₁₀⁰,      # k₋ᵒf,x₋ᵒf
    k₂₀⁰,x₂₀⁰,      # k₋ᵒs,x₋ᵒs
    k₁₂⁰,x₁₂⁰,      # kᵒa,xᵒa,
    k₂₁⁰,x₂₁⁰ = p   # k₋ᵒa,x₋ᵒa

    # @show cons[1]

    kf = bell_diss_kT(k₁₀⁰, x₁₀⁰, f) # fast diss form L-R
    ks = bell_diss_kT(k₂₀⁰, x₂₀⁰, f) # slow diss from L-R*
    ka = bell_diss_kT(k₁₂⁰, x₁₂⁰, f)   # activation from L-R to L-R*
    k_a = bell_ass_kT( k₂₁⁰, x₂₁⁰, f)  # deactivation from L-R* to L-R
    # k₋₂ = bell_diss_kT(k₋₂⁰, x₋₂⁰, f) # X dissociation from L-R-X


    uR,uA = initial_u(p,M)
    # long expression obtained from Matlab symbolics, wonder if you could get more simplified form
    df = (exp(-(t*(k_a + ka + kf + ks - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_a*uA + ka*uA + kf*uA + 2*ka*uR - ks*uA))/(2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) + (exp(-(t*(k_a + ka + kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*uA - ka*uA - kf*uA - 2*ka*uR + ks*uA))/(2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) - (exp(-(t*(k_a + ka + kf + ks - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*((k_a/2 + ka/2 + kf/2 + ks/2 - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)/ka - (k_a + ks)/ka)*(k_a/2 + ka/2 + kf/2 + ks/2 - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_a*uA + ka*uA + kf*uA + 2*ka*uR - ks*uA))/(2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) - (exp(-(t*(k_a + ka + kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*((k_a/2 + ka/2 + kf/2 + ks/2 + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)/ka - (k_a + ks)/ka)*(k_a/2 + ka/2 + kf/2 + ks/2 + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*uA - ka*uA - kf*uA - 2*ka*uR + ks*uA))/(2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))

end

function trimolecular_catch_ana_matlab(p,M,t;f=0.0)
    k₋ᵒf,x₋ᵒf,
        k₋ᵒs,x₋ᵒs,
        kᵒa,xᵒa,
        k₋ᵒa,x₋ᵒa = M.cons
    k₋₂⁰,x₋₂⁰ = p 

    kf = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
    ks = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
    ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
    k_a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)  # deactivation from L-R* to L-R
    k_2 = bell_diss_kT(k₋₂⁰, x₋₂⁰, f) # X dissociation from L-R-X
    # model

    uR,uA,uRX = initial_u(p,M)
    # long expression obtained from Matlab symbolics, wonder if you could get more simplified form
    # to weak state
    # pdf = (exp(-(t*(k_a + ka + kf + ks - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(k_2*ka^2 - 2*k_2^2*ka - k_2*k_a^2*uA + k_2^2*k_a*uA + k_2^2*ka*uA + k_2*kf^2*uA - k_2^2*kf*uA - k_a*kf^2*uA + k_a^2*kf*uA + k_2*ka^2*uR - k_2*ks^2*uA + k_2^2*ks*uA + ka*ks^2*uA - ka^2*ks*uA + kf*ks^2*uA - kf^2*ks*uA - 2*ka^2*ks*uR + k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka + k_2*ka*kf + k_2*ka*ks - k_2^2*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka*uA + k_2*ka*kf*uA - k_a*ka*kf*uA + k_2*k_a*ka*uR - 2*k_2*k_a*ks*uA - k_2*ka*ks*uA + k_a*ka*ks*uA + k_2*ka*kf*uR - 2*k_a*ka*kf*uR + 2*k_a*kf*ks*uA - 2*ka*kf*ks*uA + k_2*ka*ks*uR - 2*ka*kf*ks*uR + k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(2*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) - k_2*exp(-k_2*t)*(uA + uR - 1) + (k_2^2*exp(-k_2*t)*(k_a - k_2 + ks)*(uA + uR - 1))/(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2) + (k_2^2*ka*exp(-k_2*t)*(uA + uR - 1))/(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2) - (exp(-(t*(k_a + ka + kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(k_2*ka^2 - 2*k_2^2*ka - k_2*k_a^2*uA + k_2^2*k_a*uA + k_2^2*ka*uA + k_2*kf^2*uA - k_2^2*kf*uA - k_a*kf^2*uA + k_a^2*kf*uA + k_2*ka^2*uR - k_2*ks^2*uA + k_2^2*ks*uA + ka*ks^2*uA - ka^2*ks*uA + kf*ks^2*uA - kf^2*ks*uA - 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka + k_2*ka*kf + k_2*ka*ks + k_2^2*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka*uA + k_2*ka*kf*uA - k_a*ka*kf*uA + k_2*k_a*ka*uR - 2*k_2*k_a*ks*uA - k_2*ka*ks*uA + k_a*ka*ks*uA + k_2*ka*kf*uR - 2*k_a*ka*kf*uR + 2*k_a*kf*ks*uA - 2*ka*kf*ks*uA + k_2*ka*ks*uR - 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(2*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) + (exp(-(t*(k_a + ka + kf + ks - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(k_a - ka - kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))*(k_2*ka^2 - 2*k_2^2*ka - k_2*k_a^2*uA + k_2^2*k_a*uA + k_2^2*ka*uA + k_2*kf^2*uA - k_2^2*kf*uA - k_a*kf^2*uA + k_a^2*kf*uA + k_2*ka^2*uR - k_2*ks^2*uA + k_2^2*ks*uA + ka*ks^2*uA - ka^2*ks*uA + kf*ks^2*uA - kf^2*ks*uA - 2*ka^2*ks*uR + k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka + k_2*ka*kf + k_2*ka*ks - k_2^2*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka*uA + k_2*ka*kf*uA - k_a*ka*kf*uA + k_2*k_a*ka*uR - 2*k_2*k_a*ks*uA - k_2*ka*ks*uA + k_a*ka*ks*uA + k_2*ka*kf*uR - 2*k_a*ka*kf*uR + 2*k_a*kf*ks*uA - 2*ka*kf*ks*uA + k_2*ka*ks*uR - 2*ka*kf*ks*uR + k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(4*ka*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) + (exp(-(t*(k_a + ka + kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(ka - k_a + kf - ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))*(k_2*ka^2 - 2*k_2^2*ka - k_2*k_a^2*uA + k_2^2*k_a*uA + k_2^2*ka*uA + k_2*kf^2*uA - k_2^2*kf*uA - k_a*kf^2*uA + k_a^2*kf*uA + k_2*ka^2*uR - k_2*ks^2*uA + k_2^2*ks*uA + ka*ks^2*uA - ka^2*ks*uA + kf*ks^2*uA - kf^2*ks*uA - 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka + k_2*ka*kf + k_2*ka*ks + k_2^2*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka*uA + k_2*ka*kf*uA - k_a*ka*kf*uA + k_2*k_a*ka*uR - 2*k_2*k_a*ks*uA - k_2*ka*ks*uA + k_a*ka*ks*uA + k_2*ka*kf*uR - 2*k_a*ka*kf*uR + 2*k_a*kf*ks*uA - 2*ka*kf*ks*uA + k_2*ka*ks*uR - 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(4*ka*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))
    # to strong state
    pdf = (k_2^2*exp(-k_2*t)*(ka - k_2 + kf)*(uA + uR - 1))/(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2) - (exp(-(t*(k_a + ka + kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(k_2^2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2^2*k_a + k_2*ka^2 - k_2^2*ka + k_2*kf^2 - k_2^2*kf + k_2^2*ks - k_2*k_a^2*uA - k_a*kf^2*uA + k_a^2*kf*uA - k_2^2*k_a*uR + k_2*ka^2*uR - k_2^2*ka*uR - k_2*ks^2*uA + ka*ks^2*uA - ka^2*ks*uA - k_2*kf^2*uR + k_2^2*kf*uR + kf*ks^2*uA - kf^2*ks*uA - k_2^2*ks*uR - 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka - k_2*k_a*kf + 2*k_2*ka*kf - k_2*ka*ks - k_2*kf*ks - k_2^2*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka*uA + k_2*k_a*kf*uA - k_a*ka*kf*uA + k_2*k_a*ka*uR - 2*k_2*k_a*ks*uA + k_2*ka*ks*uA + k_a*ka*ks*uA + k_2*k_a*kf*uR - 2*k_a*ka*kf*uR + k_2*kf*ks*uA + 2*k_a*kf*ks*uA - 2*ka*kf*ks*uA + 3*k_2*ka*ks*uR + k_2*kf*ks*uR - 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(2*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) - k_2*exp(-k_2*t)*(uA + uR - 1) - (exp(-(t*(k_a + ka + kf + ks - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(k_2^2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2^2*k_a - k_2*ka^2 + k_2^2*ka - k_2*kf^2 + k_2^2*kf - k_2^2*ks + k_2*k_a^2*uA + k_a*kf^2*uA - k_a^2*kf*uA + k_2^2*k_a*uR - k_2*ka^2*uR + k_2^2*ka*uR + k_2*ks^2*uA - ka*ks^2*uA + ka^2*ks*uA + k_2*kf^2*uR - k_2^2*kf*uR - kf*ks^2*uA + kf^2*ks*uA + k_2^2*ks*uR + 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka + k_2*k_a*kf - 2*k_2*ka*kf + k_2*ka*ks + k_2*kf*ks - k_2^2*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka*uA - k_2*k_a*kf*uA + k_a*ka*kf*uA - k_2*k_a*ka*uR + 2*k_2*k_a*ks*uA - k_2*ka*ks*uA - k_a*ka*ks*uA - k_2*k_a*kf*uR + 2*k_a*ka*kf*uR - k_2*kf*ks*uA - 2*k_a*kf*ks*uA + 2*ka*kf*ks*uA - 3*k_2*ka*ks*uR - k_2*kf*ks*uR + 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(2*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) + (k_2^2*k_a*exp(-k_2*t)*(uA + uR - 1))/(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2) + (exp(-(t*(k_a + ka + kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(ka - k_a + kf - ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))*(k_2^2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2^2*k_a + k_2*ka^2 - k_2^2*ka + k_2*kf^2 - k_2^2*kf + k_2^2*ks - k_2*k_a^2*uA - k_a*kf^2*uA + k_a^2*kf*uA - k_2^2*k_a*uR + k_2*ka^2*uR - k_2^2*ka*uR - k_2*ks^2*uA + ka*ks^2*uA - ka^2*ks*uA - k_2*kf^2*uR + k_2^2*kf*uR + kf*ks^2*uA - kf^2*ks*uA - k_2^2*ks*uR - 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka - k_2*k_a*kf + 2*k_2*ka*kf - k_2*ka*ks - k_2*kf*ks - k_2^2*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka*uA + k_2*k_a*kf*uA - k_a*ka*kf*uA + k_2*k_a*ka*uR - 2*k_2*k_a*ks*uA + k_2*ka*ks*uA + k_a*ka*ks*uA + k_2*k_a*kf*uR - 2*k_a*ka*kf*uR + k_2*kf*ks*uA + 2*k_a*kf*ks*uA - 2*ka*kf*ks*uA + 3*k_2*ka*ks*uR + k_2*kf*ks*uR - 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(4*ka*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) - (exp(-(t*(k_a + ka + kf + ks - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(k_a - ka - kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))*(k_2^2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2^2*k_a - k_2*ka^2 + k_2^2*ka - k_2*kf^2 + k_2^2*kf - k_2^2*ks + k_2*k_a^2*uA + k_a*kf^2*uA - k_a^2*kf*uA + k_2^2*k_a*uR - k_2*ka^2*uR + k_2^2*ka*uR + k_2*ks^2*uA - ka*ks^2*uA + ka^2*ks*uA + k_2*kf^2*uR - k_2^2*kf*uR - kf*ks^2*uA + kf^2*ks*uA + k_2^2*ks*uR + 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka + k_2*k_a*kf - 2*k_2*ka*kf + k_2*ka*ks + k_2*kf*ks - k_2^2*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka*uA - k_2*k_a*kf*uA + k_a*ka*kf*uA - k_2*k_a*ka*uR + 2*k_2*k_a*ks*uA - k_2*ka*ks*uA - k_a*ka*ks*uA - k_2*k_a*kf*uR + 2*k_a*ka*kf*uR - k_2*kf*ks*uA - 2*k_a*kf*ks*uA + 2*ka*kf*ks*uA - 3*k_2*ka*ks*uR - k_2*kf*ks*uR + 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(4*ka*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))

end

# (k_2^2*exp(-k_2*t)*(ka - k_2 + kf)*(uA + uR - 1))/(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2) - (exp(-(t*(k_a + ka + kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(k_2^2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2^2*k_a + k_2*ka^2 - k_2^2*ka + k_2*kf^2 - k_2^2*kf + k_2^2*ks - k_2*k_a^2*uA - k_a*kf^2*uA + k_a^2*kf*uA - k_2^2*k_a*uR + k_2*ka^2*uR - k_2^2*ka*uR - k_2*ks^2*uA + ka*ks^2*uA - ka^2*ks*uA - k_2*kf^2*uR + k_2^2*kf*uR + kf*ks^2*uA - kf^2*ks*uA - k_2^2*ks*uR - 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka - k_2*k_a*kf + 2*k_2*ka*kf - k_2*ka*ks - k_2*kf*ks - k_2^2*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka*uA + k_2*k_a*kf*uA - k_a*ka*kf*uA + k_2*k_a*ka*uR - 2*k_2*k_a*ks*uA + k_2*ka*ks*uA + k_a*ka*ks*uA + k_2*k_a*kf*uR - 2*k_a*ka*kf*uR + k_2*kf*ks*uA + 2*k_a*kf*ks*uA - 2*ka*kf*ks*uA + 3*k_2*ka*ks*uR + k_2*kf*ks*uR - 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(2*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) - k_2*exp(-k_2*t)*(uA + uR - 1) - (exp(-(t*(k_a + ka + kf + ks - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(k_2^2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2^2*k_a - k_2*ka^2 + k_2^2*ka - k_2*kf^2 + k_2^2*kf - k_2^2*ks + k_2*k_a^2*uA + k_a*kf^2*uA - k_a^2*kf*uA + k_2^2*k_a*uR - k_2*ka^2*uR + k_2^2*ka*uR + k_2*ks^2*uA - ka*ks^2*uA + ka^2*ks*uA + k_2*kf^2*uR - k_2^2*kf*uR - kf*ks^2*uA + kf^2*ks*uA + k_2^2*ks*uR + 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka + k_2*k_a*kf - 2*k_2*ka*kf + k_2*ka*ks + k_2*kf*ks - k_2^2*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka*uA - k_2*k_a*kf*uA + k_a*ka*kf*uA - k_2*k_a*ka*uR + 2*k_2*k_a*ks*uA - k_2*ka*ks*uA - k_a*ka*ks*uA - k_2*k_a*kf*uR + 2*k_a*ka*kf*uR - k_2*kf*ks*uA - 2*k_a*kf*ks*uA + 2*ka*kf*ks*uA - 3*k_2*ka*ks*uR - k_2*kf*ks*uR + 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(2*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) + (k_2^2*k_a*exp(-k_2*t)*(uA + uR - 1))/(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2) + (exp(-(t*(k_a + ka + kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(ka - k_a + kf - ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))*(k_2^2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2^2*k_a + k_2*ka^2 - k_2^2*ka + k_2*kf^2 - k_2^2*kf + k_2^2*ks - k_2*k_a^2*uA - k_a*kf^2*uA + k_a^2*kf*uA - k_2^2*k_a*uR + k_2*ka^2*uR - k_2^2*ka*uR - k_2*ks^2*uA + ka*ks^2*uA - ka^2*ks*uA - k_2*kf^2*uR + k_2^2*kf*uR + kf*ks^2*uA - kf^2*ks*uA - k_2^2*ks*uR - 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka - k_2*k_a*kf + 2*k_2*ka*kf - k_2*ka*ks - k_2*kf*ks - k_2^2*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka*uA + k_2*k_a*kf*uA - k_a*ka*kf*uA + k_2*k_a*ka*uR - 2*k_2*k_a*ks*uA + k_2*ka*ks*uA + k_a*ka*ks*uA + k_2*k_a*kf*uR - 2*k_a*ka*kf*uR + k_2*kf*ks*uA + 2*k_a*kf*ks*uA - 2*ka*kf*ks*uA + 3*k_2*ka*ks*uR + k_2*kf*ks*uR - 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(4*ka*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)) - (exp(-(t*(k_a + ka + kf + ks - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/2)*(k_a/2 + ka/2 + kf/2 + ks/2 - (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)/2)*(k_a - ka - kf + ks + (k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))*(k_2^2*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2^2*k_a - k_2*ka^2 + k_2^2*ka - k_2*kf^2 + k_2^2*kf - k_2^2*ks + k_2*k_a^2*uA + k_a*kf^2*uA - k_a^2*kf*uA + k_2^2*k_a*uR - k_2*ka^2*uR + k_2^2*ka*uR + k_2*ks^2*uA - ka*ks^2*uA + ka^2*ks*uA + k_2*kf^2*uR - k_2^2*kf*uR - kf*ks^2*uA + kf^2*ks*uA + k_2^2*ks*uR + 2*ka^2*ks*uR - k_2*ka*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*kf*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*k_a*ka + k_2*k_a*kf - 2*k_2*ka*kf + k_2*ka*ks + k_2*kf*ks - k_2^2*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*k_a*ka*uA - k_2*k_a*kf*uA + k_a*ka*kf*uA - k_2*k_a*ka*uR + 2*k_2*k_a*ks*uA - k_2*ka*ks*uA - k_a*ka*ks*uA - k_2*k_a*kf*uR + 2*k_a*ka*kf*uR - k_2*kf*ks*uA - 2*k_a*kf*ks*uA + 2*ka*kf*ks*uA - 3*k_2*ka*ks*uR - k_2*kf*ks*uR + 2*ka*kf*ks*uR - k_2*k_a*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_a*kf*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*ka*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) - k_2*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + ka*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + k_2*kf*uR*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2) + kf*ks*uA*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2)))/(4*ka*(k_2*k_a + k_2*ka + k_2*kf - k_a*kf + k_2*ks - ka*ks - kf*ks - k_2^2)*(k_a^2 + 2*k_a*ka - 2*k_a*kf + 2*k_a*ks + ka^2 + 2*ka*kf - 2*ka*ks + kf^2 - 2*kf*ks + ks^2)^(1/2))


Λ₁(k₋₁,k₂,k₋₂,mₓ) = 1/2*(k₋₁ + k₋₂ + mₓ*k₂) + 1/2*√( (k₋₁ + k₋₂ + mₓ*k₂)^2 - 4*k₋₁*k₋₂ )
Λ₂(k₋₁,k₂,k₋₂,mₓ) = 1/2*(k₋₁ + k₋₂ + mₓ*k₂) - 1/2*√( (k₋₁ + k₋₂ + mₓ*k₂)^2 - 4*k₋₁*k₋₂ )

function trimolecular_ana_zhu(p,t,dens,cons)
    k₁,k₋₁,k₂,k₋₂ = p
    mᵣ,mₗ,mₓ,Ac = dens
    KₐT = k₁/k₋₁
    KₐX = k₂/k₋₂

    λ₁ = Λ₁(k₋₁,k₂,k₋₂,mₓ)
    λ₂ = Λ₂(k₋₁,k₂,k₋₂,mₓ)

    numerator = λ₂*(λ₂ - k₋₁)*exp(-λ₁*t) - λ₁*(λ₁ - k₋₁)*exp(-λ₂*t) 
    denominator = (λ₁ - λ₂)*(mₓ*k₂ + k₋₂)

    n_norm = Ac*KₐT*(1 + mₓ*KₐX)*(1 + numerator/denominator)
    # n = mᵣ*mₗ*n_norm
end

function trimolecular_ana_zhu_AcKa(p,t,dens,cons)
    AcKₐT,k₋₁,KₐX,k₋₂ = p
    mᵣ,mₗ,mₓ,Ac = dens
    # KₐT = k₁/k₋₁
    # KₐX = k₂/k₋₂
    k₂ = KₐX*k₋₂
    λ₁ = Λ₁(k₋₁,k₂,k₋₂,mₓ)
    λ₂ = Λ₂(k₋₁,k₂,k₋₂,mₓ)

    numerator = λ₂*(λ₂ - k₋₁)*exp(-λ₁*t) - λ₁*(λ₁ - k₋₁)*exp(-λ₂*t) 
    denominator = (λ₁ - λ₂)*(mₓ*k₂ + k₋₂)

    n_norm = AcKₐT*(1 + mₓ*KₐX)*(1 + numerator/denominator)
    # n = mᵣ*mₗ*n_norm
end

function trimolecular_ana_zhu_cons_bi(p,t,dens,cons)
    KₐX, k₋₂ = p
    AcKₐT, k₋₁ = cons
    mᵣ,mₗ,mₓ,Ac = dens

    k₂ = KₐX*k₋₂
    λ₁ = Λ₁(k₋₁,k₂,k₋₂,mₓ)
    λ₂ = Λ₂(k₋₁,k₂,k₋₂,mₓ)

    numerator = λ₂*(λ₂ - k₋₁)*exp(-λ₁*t) - λ₁*(λ₁ - k₋₁)*exp(-λ₂*t) 
    denominator = (λ₁ - λ₂)*(mₓ*k₂ + k₋₂)

    n_norm = AcKₐT*(1 + mₓ*KₐX)*(1 + numerator/denominator)
    # n = mᵣ*mₗ*n_norm
end

function trimolecular_ana_zhu_cons_k₋₁(p,t,dens,cons)
    AcKₐT, KₐX, k₋₂ = p
    k₋₁ = cons[1]
    mᵣ,mₗ,mₓ,Ac = dens

    k₂ = KₐX*k₋₂
    λ₁ = Λ₁(k₋₁,k₂,k₋₂,mₓ)
    λ₂ = Λ₂(k₋₁,k₂,k₋₂,mₓ)

    numerator = λ₂*(λ₂ - k₋₁)*exp(-λ₁*t) - λ₁*(λ₁ - k₋₁)*exp(-λ₂*t) 
    denominator = (λ₁ - λ₂)*(mₓ*k₂ + k₋₂)

    n_norm = AcKₐT*(1 + mₓ*KₐX)*(1 + numerator/denominator)
    # n = mᵣ*mₗ*n_norm
end

function trimolecular_ana_zhu_cons_k₋₁_k₋₂(p,t,dens,cons)
    AcKₐT, KₐX = p
    k₋₁, k₋₂ = cons
    mᵣ,mₗ,mₓ,Ac = dens

    k₂ = KₐX*k₋₂
    λ₁ = Λ₁(k₋₁,k₂,k₋₂,mₓ)
    λ₂ = Λ₂(k₋₁,k₂,k₋₂,mₓ)

    numerator = λ₂*(λ₂ - k₋₁)*exp(-λ₁*t) - λ₁*(λ₁ - k₋₁)*exp(-λ₂*t) 
    denominator = (λ₁ - λ₂)*(mₓ*k₂ + k₋₂)

    n_norm = AcKₐT*(1 + mₓ*KₐX)*(1 + numerator/denominator)
    # n = mᵣ*mₗ*n_norm
end