function cd3_dissociation!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,nₓ,n₃ = u
    k₋₃ = p[1]
    k₋₁,k₋₂ = cons
    # mᵣ,mₗ,mₓ,A = dens

    # model
    du[1] = dnᵣ = - k₋₁*nᵣ
    du[2] = dnₓ = - k₋₂*nₓ
    du[3] = dn₃ = - k₋₃*n₃

end

function bimolecular_diss_bell_activated_state_w_cons!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ,Pₐ = u
    x₋ᵒf,
        k₋ᵒs,x₋ᵒs,
        kᵒa,xᵒa,
        k₋ᵒa,x₋ᵒa = p
    k₋ᵒf = cons[1] # zero force diss, generally obtained from thermal fluct.

    k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
    ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
    k₋a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)  # deactivation from L-R* to L-R

    # model
    du[1] = dPᵣ = - k₋f.*Pᵣ - ka.*Pᵣ + k₋a.*Pₐ
    du[2] = dPₐ = + ka.*Pᵣ - k₋a.*Pₐ - k₋s.*Pₐ
end

function cd3_tri_dissociation_bell!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ₁,Pₐ₁,
        Pᵣ₂,Pₐ₂,
        Pᵣ₃,Pₐ₃ = u # 6 states

    x₋₃ᵒf,
        k₋₃ᵒs, x₋₃ᵒs,
        k₃ᵒa, x₃ᵒa,
        k₋₃ᵒa, x₋₃ᵒa = p # 7 fitted rates

    k₋₁ᵒf, x₋₁ᵒf,
        k₋₁ᵒs, x₋₁ᵒs,
        k₁ᵒa, x₁ᵒa,
        k₋₁ᵒa, x₋₁ᵒa,
        k₋₂ᵒf, x₋₂ᵒf,
        k₋₂ᵒs, x₋₂ᵒs,
        k₂ᵒa, x₂ᵒa,
        k₋₂ᵒa, x₋₂ᵒa,
        k₋₃ᵒf         = cons

    k₋₁f = bell_diss_kT(k₋₁ᵒf, x₋₁ᵒf, f) # fast diss form L-R
    k₋₁s = bell_diss_kT(k₋₁ᵒs, x₋₁ᵒs, f) # slow diss from L-R*
    k₁a  = bell_diss_kT(k₁ᵒa, x₁ᵒa, f)   # activation from L-R to L-R*
    k₋₁a = bell_ass_kT(k₋₁ᵒa, x₋₁ᵒa, f)  # deactivation from L-R* to L-R

    k₋₂f = bell_diss_kT(k₋₂ᵒf, x₋₂ᵒf, f) # fast diss form L-R
    k₋₂s = bell_diss_kT(k₋₂ᵒs, x₋₂ᵒs, f) # slow diss from L-R*
    k₂a  = bell_diss_kT(k₂ᵒa, x₂ᵒa, f)   # activation from L-R to L-R*
    k₋₂a = bell_ass_kT(k₋₂ᵒa, x₋₂ᵒa, f)  # deactivation from L-R* to L-R

    k₋₃f = bell_diss_kT(k₋₃ᵒf, x₋₃ᵒf, f) # fast diss form L-R
    k₋₃s = bell_diss_kT(k₋₃ᵒs, x₋₃ᵒs, f) # slow diss from L-R*
    k₃a  = bell_diss_kT(k₃ᵒa, x₃ᵒa, f)   # activation from L-R to L-R*
    k₋₃a = bell_ass_kT(k₋₃ᵒa, x₋₃ᵒa, f)  # deactivation from L-R* to L-R

    # model
    # γϵ
    du[1] = dPᵣ₁ = - k₋₁f*Pᵣ₁ - k₁a*Pᵣ₁ + k₋₁a*Pₐ₁
    du[2] = dPₐ₁ = + k₁a*Pᵣ₁ - k₋₁a*Pₐ₁ - k₋₁s*Pₐ₁
    # δϵ
    du[3] = dPᵣ₂ = - k₋₂f*Pᵣ₂ - k₂a*Pᵣ₂ + k₋₂a*Pₐ₂
    du[4] = dPₐ₂ = + k₂a*Pᵣ₂ - k₋₂a*Pₐ₂ - k₋₂s*Pₐ₂
    # tri
    du[5] = dPᵣ₃ = - k₋₃f*Pᵣ₃ - k₃a*Pᵣ₃ + k₋₃a*Pₐ₃
    du[6] = dPₐ₃ = + k₃a*Pᵣ₃ - k₋₃a*Pₐ₃ - k₋₃s*Pₐ₃

end

function cd3_tri_dissociation_bell_v2!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ₁,Pₐ₁,
        Pᵣ₂,Pₐ₂,
        Pᵣ₃,Pₐ₃ = u # 6 states

    # x₋₃ᵒf,
    #     k₋₃ᵒs, x₋₃ᵒs,
    #     k₃ᵒa, x₃ᵒa,
    #     k₋₃ᵒa, x₋₃ᵒa = p # 7 fitted rates
    #
    # k₋₁ᵒf, x₋₁ᵒf,
    #     k₋₁ᵒs, x₋₁ᵒs,
    #     k₁ᵒa, x₁ᵒa,
    #     k₋₁ᵒa, x₋₁ᵒa,
    #     k₋₂ᵒf, x₋₂ᵒf,
    #     k₋₂ᵒs, x₋₂ᵒs,
    #     k₂ᵒa, x₂ᵒa,
    #     k₋₂ᵒa, x₋₂ᵒa,
    #     k₋₃ᵒf         = cons

    rates = [ cons[1:2:17]; p[[2,4,6]] ]
    xs =    [ cons[2:2:16]; p[[1,3,5,7]] ]
    signs = [1,1,1,-1, 1,1,1,-1, 1,1,1,-1] # 1 is diss (+F), -1 is ass (-F)

    k₋₁f, k₋₁s, k₁a, k₋₁a,
    k₋₂f, k₋₂s, k₂a, k₋₂a,
    k₋₃f, k₋₃s, k₃a, k₋₃a = bell_diss_kT.(rates,xs,signs*f)

    # model
    # γϵ
    du[1] = dPᵣ₁ = - k₋₁f*Pᵣ₁ - k₁a*Pᵣ₁ + k₋₁a*Pₐ₁
    du[2] = dPₐ₁ = + k₁a*Pᵣ₁ - k₋₁a*Pₐ₁ - k₋₁s*Pₐ₁
    # δϵ
    du[3] = dPᵣ₂ = - k₋₂f*Pᵣ₂ - k₂a*Pᵣ₂ + k₋₂a*Pₐ₂
    du[4] = dPₐ₂ = + k₂a*Pᵣ₂ - k₋₂a*Pₐ₂ - k₋₂s*Pₐ₂
    # tri
    du[5] = dPᵣ₃ = - k₋₃f*Pᵣ₃ - k₃a*Pᵣ₃ + k₋₃a*Pₐ₃
    du[6] = dPₐ₃ = + k₃a*Pᵣ₃ - k₋₃a*Pₐ₃ - k₋₃s*Pₐ₃

end

function cd3_mix_diss_not_indep!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ₁,Pₐ₁,
        Pᵣ₂,Pₐ₂,
        Pᵣ₃ = u

    #     k₋₃₁ᵒ, x₋₃₁ᵒ,
    #     k₁₃ᵒ, x₁₃ᵒ,
    #     k₋₃₂ᵒ, x₋₃₂ᵒ,
    #     k₂₃ᵒ, x₂₃ᵒ = p

    # k₋₁ᵒf, x₋₁ᵒf,
    #     k₋₁ᵒs, x₋₁ᵒs,
    #     k₁ᵒa, x₁ᵒa,
    #     k₋₁ᵒa, x₋₁ᵒa,
    #     k₋₂ᵒf, x₋₂ᵒf,
    #     k₋₂ᵒs, x₋₂ᵒs,
    #     k₂ᵒa, x₂ᵒa,
    #     k₋₂ᵒa, x₋₂ᵒa = cons

    rates = [ cons[1:2:15]; p[1:2:7] ]
    xs =    [ cons[2:2:16]; p[2:2:8] ]
    signs = [1,1,1,-1, 1,1,1,-1, 1,-1,1,-1] # 1 is diss (+F), -1 is ass (-F)

    k₋₁f, k₋₁s, k₁a, k₋₁a,
    k₋₂f, k₋₂s, k₂a, k₋₂a,
    k₋₃₁, k₁₃, k₋₃₂, k₂₃ = bell_diss_kT.(rates,xs,signs*f)

    # model
    # γϵ
    du[1] = dPᵣ₁ = - k₋₁f*Pᵣ₁ - k₁a*Pᵣ₁ + k₋₁a*Pₐ₁ - k₁₃*Pᵣ₁ + k₋₃₁*Pᵣ₃
    du[2] = dPₐ₁ = + k₁a*Pᵣ₁ - k₋₁a*Pₐ₁ - k₋₁s*Pₐ₁
    # δϵ
    du[3] = dPᵣ₂ = - k₋₂f*Pᵣ₂ - k₂a*Pᵣ₂ + k₋₂a*Pₐ₂ - k₂₃*Pᵣ₂ + k₋₃₂*Pᵣ₃
    du[4] = dPₐ₂ = + k₂a*Pᵣ₂ - k₋₂a*Pₐ₂ - k₋₂s*Pₐ₂
    # tri
    du[5] = dPᵣ₃ = - k₋₃₁*Pᵣ₃ - k₋₃₂*Pᵣ₃ + k₁₃*Pᵣ₁ + k₂₃*Pᵣ₂

end

function tri_diss_bell_r!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ,Pₐ,P₃ = u
    k₋ᵒf,x₋ᵒf,
        k₋ᵒs,x₋ᵒs,
        kᵒa,xᵒa,
        k₋ᵒa,x₋ᵒa = cons
    k₂⁰,x₂⁰,  
        k₋₂⁰,x₋₂⁰ = p 

    k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
    ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
    k₋a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)  # deactivation from L-R* to L-R
    k₂  = bell_ass_kT(k₂⁰, x₂⁰, f)   # X association w/ L-R
    k₋₂ = bell_diss_kT(k₋₂⁰, x₋₂⁰, f) # X dissociation from L-R-X
    # model
    du[1] = dPᵣ = - k₋f.*Pᵣ - ka.*Pᵣ + k₋a.*Pₐ - k₂.*Pᵣ + k₋₂.*P₃
    du[2] = dPₐ = + ka.*Pᵣ - k₋a.*Pₐ - k₋s.*Pₐ
    du[3] = dP₃ = + k₂.*Pᵣ - k₋₂.*P₃
end

function tri_diss_bell_a!(du,u,p,t,dens,cons;f=0.0)
    # unpack rates and constants
    Pᵣ,Pₐ,P₃ = u
    k₋ᵒf,x₋ᵒf,
        k₋ᵒs,x₋ᵒs,
        kᵒa,xᵒa,
        k₋ᵒa,x₋ᵒa = cons
    k₂⁰,x₂⁰,  
        k₋₂⁰,x₋₂⁰ = p 

    k₋f = bell_diss_kT(k₋ᵒf, x₋ᵒf, f) # fast diss form L-R
    k₋s = bell_diss_kT(k₋ᵒs, x₋ᵒs, f) # slow diss from L-R*
    ka  = bell_diss_kT(kᵒa, xᵒa, f)   # activation from L-R to L-R*
    k₋a = bell_ass_kT(k₋ᵒa, x₋ᵒa, f)  # deactivation from L-R* to L-R
    k₂  = bell_ass_kT(k₂⁰, x₂⁰, f)   # X association w/ L-R
    k₋₂ = bell_diss_kT(k₋₂⁰, x₋₂⁰, f) # X dissociation from L-R-X
    # model
    du[1] = dPᵣ = - k₋f.*Pᵣ - ka.*Pᵣ + k₋a.*Pₐ 
    du[2] = dPₐ = + ka.*Pᵣ - k₋a.*Pₐ - k₋s.*Pₐ - k₂.*Pₐ + k₋₂.*P₃
    du[3] = dP₃ = + k₂.*Pₐ - k₋₂.*P₃
end