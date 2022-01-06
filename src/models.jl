function bimolecular!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ = u[1]
    k₁,k₋₁  = p
    mᵣ,mₗ,A = dens
    # model
    du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ

end


function trimolecular!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,n₃ = u
    k₂,k₋₂ = p
    k₁,k₋₁ = cons
    # k₁,k₋₁ = [1e-4,0.5]
    mᵣ,mₗ,mₓ,A = dens
    # model
    du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ - A*k₂*mₓ*nᵣ + k₋₂*n₃
    du[2] = dn₃ = A*k₂*mₓ*nᵣ - k₋₂*n₃

end

function trimolecular_no_cons!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,n₃ = u
    k₁,k₋₁,k₂,k₋₂ = p
    # k₁,k₋₁ = cons
    # k₁,k₋₁ = [1e-4,0.5]
    mᵣ,mₗ,mₓ,A = dens
    # model
    du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ - A*k₂*mₓ*nᵣ + k₋₂*n₃
    du[2] = dn₃ = A*k₂*mₓ*nᵣ - k₋₂*n₃

end

function two_path!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,nₓ,n₃ = u
    k₃ ,k₋₃, k₄ ,k₋₄ = p
    k₁,k₋₁,k₂,k₋₂ = cons
    mᵣ,mₗ,mₓ,A = dens

    # model
    du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ - A*k₄*mₓ*nᵣ + k₋₄*n₃
    du[2] = dnₓ = A*k₂*mₓ*mₗ - k₋₂*nₓ - A*k₃*mᵣ*nₓ + k₋₃*n₃
    du[3] = dn₃ = A*k₄*mₓ*nᵣ + A*k₃*mᵣ*nₓ - k₋₄*n₃ - k₋₃*n₃

end

const kT = 4.114 # pN*nm

bell_diss(k₋ᵒ,fᵒ,f) = k₋ᵒ .* exp(f/fᵒ)
bell_ass(k₊ᵒ,fᵒ,f) = k₊ᵒ .* exp(-f/fᵒ)

bell_diss_kT(k₋ᵒ,x,f) = k₋ᵒ .* exp(f .* x ./ kT)
bell_ass_kT(k₊ᵒ,x,f) = k₊ᵒ .* exp(-f .*x ./ kT)

function bimolecular_diss!(du,u,p,t,dens,cons)
    # unpack rates and constants
    Pᵣ = u[1]
    k₋₁  = p[1]
    # model
    du[1] = dPᵣ = - k₋₁.*Pᵣ

end

function trimolecular_diss!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,n₃ = u
    k₂,k₋₂ = p
    # k₂ = 0.0
    k₋₁ = cons[1]
    mᵣ,mₗ,mₓ,A = dens
    # model
    du[1] = dnᵣ = - k₋₁*nᵣ - A*k₂*mₓ*nᵣ + k₋₂*n₃
    du[2] = dn₃ = A*k₂*mₓ*nᵣ - k₋₂*n₃

end

function trimolecular_diss_no_on2!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,n₃ = u
    k₋₂ = p[1]
    # k₂ = 0.0
    k₋₁ = cons[1]
    mᵣ,mₗ,mₓ,A = dens
    # model
    du[1] = dnᵣ = - k₋₁*nᵣ + k₋₂*n₃
    du[2] = dn₃ = - k₋₂*n₃

end

function bimolecular_diss_two_states!(du,u,p,t,dens,cons)
    # unpack rates and constants
    # P = u[1]
    Pᵣ,Pₐ = u
    k₋f,k₋s = p
    # k₋f = cons[1]

    # model
    du[1] = dPᵣ = - k₋f.*Pᵣ
    du[2] = dPₐ = - k₋s.*Pₐ
end

function bimolecular_diss_bell!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ,Pₐ = u
    k₋ᵒf,x₋ᵒf,
        k₋ᵒs,x₋ᵒs,
        kᵒa,xᵒa,
        k₋ᵒa,x₋ᵒa = p
    # f = cons[1]

    k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
    ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
    k₋a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)  # deactivation from L-R* to L-R

    # model
    du[1] = dPᵣ = - k₋f.*Pᵣ - ka.*Pᵣ + k₋a.*Pₐ
    du[2] = dPₐ = + ka.*Pᵣ - k₋a.*Pₐ - k₋s.*Pₐ
end


function bimolecular_diss_bell_single!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ = u[1]
    k₋ᵒf,x₋ᵒf = p
    # f = cons[1]

    k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R

    # model
    du[1] = dPᵣ = - k₋f.*Pᵣ
end

function bimolecular_diss_bell_w_cons!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ = u[1]
    x₋ᵒ = p[1]

    k₋ᵒ = cons[1]
    # f = cons[2]

    k₋ = bell_diss_kT(k₋ᵒ, x₋ᵒ, f) # fast diss form L-R

    # model
    du[1] = dPᵣ = - k₋.*Pᵣ
end

function bimolecular_diss_two_path!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    # P = u[1]
    Pᵣ,Pₐ = u
    k₋ᵒf,x₋ᵒf,k₋ᵒs,x₋ᵒs, = p

    k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*

    # model
    du[1] = dPᵣ = - k₋f.*Pᵣ
    du[2] = dPₐ = - k₋s.*Pₐ

    # du[1] = dPᵣ + dPₐ
end

function trimolecular_diss_two_path!(du,u,p,t,dens,cons)
    # unpack rates and constants
    # P = u[1]
    Pᵣ,Pₐ,P₃ = u
    k₋ᵒf,x₋ᵒf,k₋ᵒs,x₋ᵒs,kᵒ3,xᵒ3,
        k₋ᵒ3,x₋ᵒ3 = p
    f = cons[1]

    k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
    k₃ = bell_ass_kT(kᵒ3, xᵒ3, f)
    k₋₃ = bell_diss_kT(k₋ᵒ3, x₋ᵒ3, f)

    du[1] = dPᵣ = - k₋f.*Pᵣ - k₃.*Pᵣ + (1-w).*k₋₃.*P₃
    du[2] = dPₐ = - k₋s.*Pₐ - k₃.*Pₐ +    w .*k₋₃.*P₃
    du[3] = dP₃ = + k₃.*(Pᵣ+Pₐ) - k₋₃.*P₃

end

function trimolecular_diss_bell!(du,u,p,t,dens,cons)
    # unpack rates and constants
    Pᵣ,Pₐ,P₃ = u

    kᵒ3,xᵒ3,
    k₋ᵒ3,x₋ᵒ3 = p

    f = cons[1]

    k₋ᵒf,x₋ᵒf,
        k₋ᵒs,x₋ᵒs,
        kᵒa,xᵒa,
        k₋ᵒa,x₋ᵒa  = bimolecular_cons

    k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
    ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
    k₋a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)  # deactivation from L-R* to L-R

    k₃ = bell_ass_kT(kᵒ3, xᵒ3, f)
    k₋₃ = bell_diss_kT(k₋ᵒ3, x₋ᵒ3, f)


    # model
    du[1] = dPᵣ = - k₋f.*Pᵣ - ka.*Pᵣ + k₋a.*Pₐ - k₃.*Pᵣ + k₋₃.*P₃
    du[2] = dPₐ = + ka.*Pᵣ - k₋a.*Pₐ - k₋s.*Pₐ
    du[3] = dP₃ = + k₃.*Pᵣ - k₋₃.*P₃
end


function bimolecular_catch_two_states!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ,Pₐ = u
    k₋ᵒc,x₋ᵒc,
        k₋ᵒs,x₋ᵒs = p

    k₋c = bell_ass_kT(k₋ᵒc, x₋ᵒc, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*

    # model
    du[1] = dPᵣ = - k₋c.*Pᵣ
    du[2] = dPₐ = - k₋s.*Pₐ
end

function bimolecular_catch_1_state_2_paths!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ = u[1]
    k₋ᵒc,x₋ᵒc,
        k₋ᵒs,x₋ᵒs = p

    k₋c = bell_ass_kT(k₋ᵒc, x₋ᵒc, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*

    # model
    du[1] = dPᵣ = - k₋c.*Pᵣ - k₋s.*Pᵣ
end


## functions to calculate initial bond distribution vector (at t=0)

#=
if the delay between binding and beginning of the measurement is sufficiently
long, the probability of occupying either state will reach an equilibrium. This
equilibrium is independent of force
=#
function initial_state(kᵒa,k₋ᵒa)
        u₀1 = k₋ᵒa / (kᵒa + k₋ᵒa)
        u₀ = [u₀1, (1.0 - u₀1)]
end

#=
assumes that the time between binding and measurement is short relative to the
transition rates between the two bound states, then the bond will stay in the
initial state it occupied upon binding. The probabilities of occupying each of
these states can be derived from the principle of detailed balance.
=#
function initial_state(k₋ᵒf,k₋ᵒs,kᵒa,k₋ᵒa)
        u₀1 = (k₋ᵒa * k₋ᵒf) / (k₋ᵒa * k₋ᵒf + kᵒa * k₋ᵒs)
        u₀ = [u₀1, (1.0 - u₀1)]
end

#=
assumes enough time to reach equilibrium, and that the bond experiences some
fraction a of the measured force.
IMPORTANT: must add a*f to function call: initial_state(p, a * f)
=#


function initial_state(p::Vector,f)
        k₋ᵒf,x₋ᵒf,
        k₋ᵒs,x₋ᵒs,
        kᵒa,xᵒa,
        k₋ᵒa,x₋ᵒa = p

        k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
        k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
        ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
        k₋a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)

        u₀1 = k₋a / (ka + k₋a)
        u₀ = [u₀1, (1.0 - u₀1)]
end

# function
# u = []
# for f in Fs
#         u₀ = initial_state(rates,f)
#         push!(u,u₀)
# end
#
# plot(reduce(hcat,u))

function initial_u(k,M)
    if ~isempty(M.cons)
        k = vcat(M.cons,k)
    end

    if M.u₀_type == "equilibrium"
        k₋ᵒf,x₋ᵒf,k₋ᵒs,x₋ᵒs,kᵒa,xᵒa,k₋ᵒa,x₋ᵒa = k[end-7:end]
        new_u₀ = initial_state(kᵒa,k₋ᵒa)

    elseif M.u₀_type == "detailed_balance"
        k₋ᵒf,x₋ᵒf,k₋ᵒs,x₋ᵒs,kᵒa,xᵒa,k₋ᵒa,x₋ᵒa = k[end-7:end]
        new_u₀ = initial_state(k₋ᵒf,k₋ᵒs,kᵒa,k₋ᵒa)

    elseif M.u₀_type == "force"
        f = M.cons
        a = k[end]
        new_u₀ = initial_state(k,a*f)

    else
        new_u₀ = []
    end

    if isempty(M.u₀)
        u₀ = new_u₀
    else
        u₀ = M.u₀
        u₀ = convert.(eltype(new_u₀),u₀)
        u₀ = append!(u₀, new_u₀) # append new initial states to provided initial states
    end
    return u₀
end

