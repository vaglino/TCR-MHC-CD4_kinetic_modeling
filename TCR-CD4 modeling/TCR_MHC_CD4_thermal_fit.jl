## CD3 γϵ and δϵ cooperative binding to TCRαβ ecotodomain

# bimolecular dissociations
#       γϵ + TCR ← γϵ-TCR
#       δϵ + TCR ← δϵ-TCR
##
using Revise
# include("model_fit_global.jl")
dir = "C:/Users/stravaglino3/Documents/TCR-MHC-CD4 models/src/"
includet(dir*"models.jl")
includet(dir*"models_dissociation.jl")
includet(dir*"helper_functions.jl")
includet(dir*"binning_helper_fs.jl")
include(dir*"fitting.jl")
includet(dir*"mle_fitting.jl")
includet(dir*"helper_plotting.jl")
using DataFrames, Plots, Measurements
##

## load thermal data
data_dir = "C:/Users/stravaglino3/Documents/TCR-MHC-CD4 models/TCR_CD4_project/dissociation data/"
thermal_b = CSV.read(data_dir*"Thermal fluctuation assay _ CD4-E8TCR-TPI_HLADR.csv",DataFrame)

t_tm_b_thrml = delete_missing(thermal_b[:,3:3])
t_cm_b_thrml = delete_missing(thermal_b[:,8:8])
t_tmc_b_thrml = delete_missing(thermal_b[:,13:13])
p_tm_b_thrml = survival_p.(t_tm_b_thrml)
p_cm_b_thrml = survival_p.(t_cm_b_thrml)
p_tmc_b_thrml = survival_p.(t_tmc_b_thrml)


tm_b_thrml = DissData(p_tm_b_thrml,t_tm_b_thrml,[[]])
cm_b_thrml = DissData(p_cm_b_thrml,t_cm_b_thrml,[[]])
tmc_b_thrml = DissData(p_tmc_b_thrml,t_tmc_b_thrml,[[]])

# x_min = minimum(t_mix_thrml[1])
# γϵ_thrml = DissData(p_γϵ_thrml,rescale_t.(t_γϵ_thrml;Δx=x_min),[[]])
# δϵ_thrml = DissData(p_δϵ_thrml,rescale_t.(t_δϵ_thrml;Δx=x_min),[[]])
# mix_thrml = DissData(p_mix_thrml,t_mix_thrml,[[]])


## Thermal fluctuation bimolecular dissociation

u₀ = [1.0]           # Initial condition
tspan = (0.0,10.0)  # Simulation interval
dens =  fill([],length(t_tm_b_thrml))
model = bimolecular_diss!
cons_thrml = [[0.0]]

# pack parameters
tm_b_M_thrml = M_global(model,u₀,tspan,dens,cons_thrml)

bounds = [(0.,100.)]
pinit = [0.0]

## CD3γϵ thermal fitting---------------------------------------------------------
res_tm_b_thrml = global_LS_fit(tm_b_M_thrml, tm_b_thrml, bounds)
k₋₁ᵒ = res_tm_b_thrml.rates[1]
@show k_tm_b_thrml = k₋₁ᵒ ±  √(res_tm_b_thrml.hes[1])
@show res_tm_b_thrml.l

## CD3δϵ thermal fitting---------------------------------------------------------
cm_b_M_thrml = M_global(model,u₀,tspan,dens,cons_thrml)
res_cm_b_thrml = global_LS_fit(cm_b_M_thrml, cm_b_thrml, bounds)
k₋₂ᵒ = res_cm_b_thrml.rates[1]
@show k_cm_b_thrml = k₋₂ᵒ ±  √(res_cm_b_thrml.hes[1])
@show res_cm_b_thrml.l

## CD3 mix thermal fitting-------------------------------------------------------
tmc_b_M_thrml_1 = M_global(model,u₀,tspan,dens,cons_thrml)
res_tmc_b_thrml_1 = global_LS_fit(tmc_b_M_thrml_1, tmc_b_thrml, bounds)
k₋₃ᵒ_1 = res_tmc_b_thrml_1.rates[1]
@show k_tmc_b_thrml_1 = k₋₃ᵒ_1 ± √(res_tmc_b_thrml_1.hes[1])
@show res_tmc_b_thrml_1.l

# Kₐγϵ, Kₐδϵ = [2.6e-6, 2.1e-6] # bimolecular affinities from adhesion frequency assay
# u₀_bi = [Kₐγϵ/(Kₐγϵ + Kₐδϵ), Kₐδϵ/(Kₐγϵ + Kₐδϵ)] # initial proportion of bimolecular bonds from relative Ka
# u₀ = u₀_bi
# cons_mix_thrml = [[k₋₁ᵒ[1],k₋₂ᵒ[1]]] # keep bimolecuar rates constants
# cons_mix_thrml = [[1/mean(γϵ_thrml.t[1]),1/mean(δϵ_thrml.t[1])]] # keep bimolecular rates constants
# cons_tmc_b_thrml = k₋₁ᵒ
# u₀ = [0.9,0.1]
# model = trimolecular_diss!
# dens_tmc_b = [[25,38,20,1]]
# bounds = [(0.0,10.0), (0.0,10.0)] # [k2,k_2]
# pinit = [1.0,1.0]
# tmc_b_M_thrml = M_global(model,u₀,tspan,dens_tmc_b,cons_tmc_b_thrml)

# res_tmc_b_thrml = global_LS_fit(tmc_b_M_thrml, tmc_b_thrml, bounds)
# k₋₃ᵒ = res_tmc_b_thrml.rates
# σ₃ = hessian2σ_LS(res_tmc_b_thrml.hes,res_tmc_b_thrml.l,length(tmc_b_thrml.t[1]),2);
# k_tmc_b_thrml = res_tmc_b_thrml.rates .± σ₃


# 2 independent dissociations
# cons_tmc_b_thrml = k₋₁ᵒ
u₀ = [0.9,0.1]
model = bimolecular_diss_two_states!
dens_tmc_b = [[25,38,20,1]]
bounds = [(0.0,10.0), (0.0,10.0), (0.0,0.4)] # [k2,k_2]
pinit = [1.0,0.1,0.3]
tmc_b_M_thrml_2 = M_global(model,u₀,tspan,dens_tmc_b,[[]])
res_tmc_b_thrml_2 = global_LS_fit(tmc_b_M_thrml_2, tmc_b_thrml, bounds)
k₋₃ᵒ_2 = res_tmc_b_thrml_2.rates
@show res_tmc_b_thrml_2.l


# sequential TCR-MHC-CD4 -> TCR-MHC + CD4 -> ⌽
cons_tmc_b_thrml = k₋₁ᵒ
u₀ = [0.9,0.1]
model = trimolecular_diss!
dens_tmc_b = [[25,38,20,1]]
bounds = [(0.0,100.0), (0.0,100.0), (0.0,1.0) ] # [k2,k_2]
pinit = [1.0,1.0,0.1]
tmc_b_M_thrml_3 = M_global(model,u₀,tspan,dens_tmc_b,cons_tmc_b_thrml)
res_tmc_b_thrml_3 = global_LS_fit(tmc_b_M_thrml_3, tmc_b_thrml, bounds)
@show k₋₃ᵒ_3 = res_tmc_b_thrml_3.rates
@show res_tmc_b_thrml_3.l
# k₋₃ᵒ_3 = [ 0, 5.045339922925327, 0.2110701531504992]

cons_tmc_b_thrml = k₋₁ᵒ
u₀ = [0.9,0.1]
model = trimolecular_diss_no_on2!
dens_tmc_b = [[25,38,20,1]]
bounds = [(0.0,100.0), (0.0,0.3) ] # [k2,k_2]
pinit = [1.0,0.1]
tmc_b_M_thrml_3 = M_global(model,u₀,tspan,dens_tmc_b,cons_tmc_b_thrml)
res_tmc_b_thrml_3 = global_LS_fit(tmc_b_M_thrml_3, tmc_b_thrml, bounds)
@show k₋₃ᵒ_3 = res_tmc_b_thrml_3.rates
@show res_tmc_b_thrml_3.l
# k₋₃ᵒ_3 = [ 5.045339922925327, 0.2110701531504992]

# sequential TCR-MHC-CD4 -> CD4-MHC + TCR -> ⌽
cons_tmc_b_thrml = k₋₂ᵒ
u₀ = [0.9,0.1]
model = trimolecular_diss!
# dens_tmc_b = [[25,38,20,1]]
dens_tmc_b = [[20,38,25,1]]

bounds = [(0.0,100.0), (0.0,100.0), (0.0,1.0) ] # [k2,k_2]
pinit = [1.0,1.0,0.1]
tmc_b_M_thrml_4 = M_global(model,u₀,tspan,dens_tmc_b,cons_tmc_b_thrml)
res_tmc_b_thrml_4 = global_LS_fit(tmc_b_M_thrml_4, tmc_b_thrml, bounds)
@show k₋₃ᵒ_4 = res_tmc_b_thrml_4.rates
@show res_tmc_b_thrml_4.l



# u₀_tri = rates_mix_thrml[1]
# k₋₃ᵒ = rates_mix_thrml[2]

#--------------------------------------------------------------------------------
## thermal fluctuation plotting
# pyplot() #backend
plotly()
# gr()
p_thrml = plot()
fit_thrml_tm_b, p1 = show_diss_fit!(p_thrml,tm_b_M_thrml,k₋₁ᵒ, tm_b_thrml; cmap=:blue,label="TCR-MHC")
fit_thrml_cm_b, p2 = show_diss_fit!(p_thrml,cm_b_M_thrml,k₋₂ᵒ, cm_b_thrml; cmap=:red,label="CD4-MHC")
fit_thrml_tmc_b_1, p3_1 = show_diss_fit!(p_thrml,tmc_b_M_thrml_1,k₋₃ᵒ_1, tmc_b_thrml; cmap=:green,label="TCR-MHC-CD4")
fit_thrml_tmc_b_2, p3_2 = show_diss_fit!(p_thrml,tmc_b_M_thrml_2,k₋₃ᵒ_2, tmc_b_thrml; cmap=:green,label="TCR-MHC-CD4")
fit_thrml_tmc_b_3, p3_3 = show_diss_fit!(p_thrml,tmc_b_M_thrml_3,k₋₃ᵒ_3, tmc_b_thrml; cmap=:green,label="TCR-MHC-CD4")
fit_thrml_tmc_b_4, p3_4 = show_diss_fit!(p_thrml,tmc_b_M_thrml_4,k₋₃ᵒ_4, tmc_b_thrml; cmap=:green,label="TCR-MHC-CD4")

plot(p1)
xlims!(0,1.2)
ylims!(-5,0.1)
tmp_sol = run_model(mixM_thrml,rates_mix_thrml,mixM_thrml.dens,mixM_thrml.cons[1],0:0.1:10) 
mean_lifetime_from_S(tmp_sol)
summed_sol = vec(sum(Array(tmp_sol),dims=1))


y_thrml_γϵ = vec(get_y_thrml(fit_thrml_γϵ[1]))
y_thrml_δϵ = vec(get_y_thrml(fit_thrml_δϵ[1]))
y_thrml_mix = vec(get_y_thrml(fit_thrml_mix[1]))

# CSV.write(results_dir*"l  np_thrml_γϵ.csv",measurements2Df(fit_thrml_γϵ[1].t,y_thrml_γϵ) )
# CSV.write(results_dir*"lnp_thrml_δϵ.csv",measurements2Df(fit_thrml_δϵ[1].t,y_thrml_δϵ) )
# CSV.write(results_dir*"lnp_thrml_mix.csv",measurements2Df(fit_thrml_mix[1].t,y_thrml_mix) )

measurements2Df.(reduce(vcat,fit_thrml_γϵ[1].u))

reduce(vcat,fit_thrml_γϵ[1].u)

function get_y_thrml(sol)
    Σ = sum(Array(sol),dims=1)
    Σ = Im2Re(Σ)
end