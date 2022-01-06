## modeling of TCR-MHC-CD4 kinetics from adhesion frequency experiments
## load dependencies
using Revise
# change main directory accordingly ro where software is saved on own machine
main_dir = "C:/Users/stravaglino3/Documents/TCR-MHC-CD4 modeling for github/"

src_dir = main_dir*"src/"
includet(src_dir*"model_fit_global_v2.jl")
includet(src_dir*"models.jl")
includet(src_dir*"analytical_models.jl")
includet(src_dir*"helper_functions.jl")
includet(src_dir*"global_plotting.jl")
using Plots
using BlackBoxOptim, Measurements, FiniteDiff

#______________________________________________________________
## Site densities (from flow cytometry)
densities_T =   [25.0, 38.0, 1.0]       # TCR, pMHC, Ac
densities_TC1 = [15.0, 38.0, 10.0, 1.0] # TCR, pMHC, CD4, Ac
densities_TC2 = [8.0,  38.0, 17.0, 1.0]
densities_TC3 = [3.0,  38.0, 22.0, 1.0]
densities_C =   [1700., 850., 1]       # CD4, pMHC, A

#______________________________________________________________
## Importing experimental data
data_dir = main_dir*"TCR-CD4 modeling/data/"

using CSV
tcr_data = CSV.read(data_dir*"tcr_only_data.csv", DataFrame)
tcr_cd4_1 = CSV.read(data_dir*"tcr_cd4_1.csv", DataFrame)
tcr_cd4_2 = CSV.read(data_dir*"tcr_cd4_2.csv", DataFrame)
tcr_cd4_3 = CSV.read(data_dir*"tcr_cd4_3.csv", DataFrame)
cd4_data = CSV.read(data_dir*"cd4_data_from_muaz.csv", DataFrame)
t_t,   Pa_t,   n_se_t   = clean_data(tcr_data)
t_tc1, Pa_tc1, n_se_tc1 = clean_data(tcr_cd4_1)
t_tc2, Pa_tc2, n_se_tc2 = clean_data(tcr_cd4_2)
t_tc3, Pa_tc3, n_se_tc3 = clean_data(tcr_cd4_3)
t_c,   Pa_c,   n_se_c   = clean_data(cd4_data)

t_t,   n_t,   n_se_t   = clean_and_normalize(tcr_data,  densities_T)
t_tc1, n_tc1, n_se_tc1 = clean_and_normalize(tcr_cd4_1, densities_TC1)
t_tc2, n_tc2, n_se_tc2 = clean_and_normalize(tcr_cd4_2, densities_TC2)
t_tc3, n_tc3, n_se_tc3 = clean_and_normalize(tcr_cd4_3, densities_TC3)
t_c,   n_c,   n_se_c   = clean_and_normalize(cd4_data,  densities_C)

# n_t   = Pa2n(Pa_t)   ./ (25.0 * 38.0); n_se_t = sem()
# n_tc1 = Pa2n(Pa_tc1) ./ (15.0 * 38.0)
# n_tc2 = Pa2n(Pa_tc2) ./ (8.0 * 38.0)
# n_tc3 = Pa2n(Pa_tc3) ./ (3.0 * 38.0)
# n_c   = Pa2n(Pa_c)   ./ (1700 * 850.0)

#______________________________________________________________
## plot experimental data
plot()
scatter(n_se_t[:,1],n_se_t[:,2],yerror=n_se_t[:,3],markercolor=:black,markersize=5)
scatter!(n_se_tc1[:,1],n_se_tc1[:,2],yerror=n_se_tc1[:,3],markercolor=:blue,markersize=5)
scatter!(n_se_tc2[:,1],n_se_tc2[:,2],yerror=n_se_tc2[:,3],markercolor=:red,markersize=5)
scatter!(n_se_tc3[:,1],n_se_tc3[:,2],yerror=n_se_tc3[:,3],markercolor=:cyan,markersize=5,legend=false)
scatter!(n_se_c[:,1],n_se_c[:,2],yerror=n_se_c[:,3],markercolor=:green,markersize=5,legend=false)

#______________________________________________________________
# First, we are going to fit both TCR-MHC and CD4-MHC bimolecular interactions for reference
## TCR-MHC fit (bimolecular)
tspan = (0,maximum(t_t)+1)
u₀ = [0.0]
pinit = [0.0, 0.0] # [k₁, k₋₁]
bounds = [(1.0e-6,100.0), (1.0e-6,100.0)] # [k₁, k₋₁]


tcr = Data([n_t],[t_t])
M_tcr = M_global(bimolecular!,u₀,tspan,[densities_T],[[]],0.0,Rodas5)

res_tcr = global_LS_fit(M_tcr,tcr,bounds)
rates_tcr = res_tcr.rates
σ_tcr = hessian2σ(res_tcr.hes,res_tcr.l,tcr.n[1],rates_tcr)
@show k_t = rates_tcr .± σ_tcr

#______________________________________________________________
## CD4-MHC fit (bimolecular)
tspan = (0,maximum(t_t)+1)
u₀ = [0.0]
bounds = [(1.0e-9,100.0),(1.0e-9,100.0)] # [k₂, k₋₂]
pinit = [0.0,0.0] # [k₂, k₋₂]

cd4 = Data([n_c],[t_c])
M_cd4 = M_global(bimolecular!,u₀,tspan,[densities_C],[[]],0.0,Rodas5)

res_cd4 = global_LS_fit(M_cd4,cd4,bounds)
rates_cd4 = res_cd4.rates
σ_cd4 = hessian2σ(res_cd4.hes,res_cd4.l,cd4.n[1],rates_cd4)
@show k_c = rates_cd4 .± σ_cd4

#______________________________________________________________
## plot bimolecular fits
plot()
scatter(n_se_t[:,1],n_se_t[:,2],yerror=n_se_t[:,3],markercolor=:black,markersize=5,label="TCR-MHC")
scatter!(n_se_c[:,1],n_se_c[:,2],yerror=n_se_c[:,3],markercolor=:green,markersize=5,label="CD4-MHC",legend=:bottomright)

plot_n_fit(M_tcr,rates_tcr,densities_T,[];Color=:black)
plot_n_fit(M_cd4,rates_cd4,densities_C,[];Color=:green)

f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_T, [])
Color = :black; fit_t = plot_sum(rates_tcr,:black)

f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_C, [])
Color = :green; fit_c = plot_sum(rates_cd4,:green)

#______________________________________________________________
## Now fit all trimolecular TCR-MHC-CD4 interactions at the same time (global fit)

densities_T_full = [25.0, 38.0, 0.0, 1.0] # TCR, pMHC, CD4, Ac
u₀ = [0.0, 0.0] # [n_TCR, n_3]
cons_tcr = []
cons_tc = fill(cons_tcr,4) # don't hold any rate constant
dens_tc = [densities_T_full,densities_TC1,densities_TC2,densities_TC3]
tspan = (0,maximum(t_tc1)+1)
# pinit = [0.0,0.0,0.0,0.0] # [k₁, k₋₁, k₂, k₋₂]
pinit = [8.25E-04, # from Dr. Zhu's excel sheet
        0.536,
        5.55E-02,
        0.975]
pinit = [1e-3, # from Dr. Zhu's excel sheet
        0.1,
        0.01,
        1]

tc = Data([n_t,n_tc1,n_tc2,n_tc3],[t_t,t_tc1,t_tc2,t_tc3]) # fit all experiments at same time

# model = trimolecular_no_cons!
model = trimolecular_ana_zhu_AcKa
M_tc = M_global(model,u₀,tspan,dens_tc,cons_tc,0.0,Rodas5)

bounds = [(1.0e-5,100.0),(1.0e-9,100.0), # [AcKaTCR, k₋₁]
        (1.0e-5, 10),(1.0e-5,2)] # [KaCD4, k₋₂]
# uncomment p₀_bbo = pinit in model_fit_global_v2
res_tc = global_LS_fit(M_tc,tc,bounds)
rates_tc = res_tc.rates
res_tc.hes
σ_tc = hessian2σ(res_tc.hes,res_tc.l,reduce(vcat,tc.n),res_tc.rates) #./ √(64)
k_tc = rates_tc .± σ_tc

#______________________________________________________________
## plotting
rates = k_tc
# rates = rates_tc
# rates = vcat(rates_tc[1:2],0.0541125,0.975)

plot()
fit_t   = plot_n_fit(M_tc,rates,densities_T_full,[];Color=:black)
fit_tc1 = plot_n_fit(M_tc,rates,densities_TC1,[];Color=:blue)
fit_tc2 = plot_n_fit(M_tc,rates,densities_TC2,[];Color=:red)
fit_tc3 = plot_n_fit(M_tc,rates,densities_TC3,[];Color=:cyan)

scatter!(n_se_t[:,1],  n_se_t[:,2],  yerror=n_se_t[:,3],  markercolor=:black,markersize=5)
scatter!(n_se_tc1[:,1],n_se_tc1[:,2],yerror=n_se_tc1[:,3],markercolor=:blue,markersize=5)
scatter!(n_se_tc2[:,1],n_se_tc2[:,2],yerror=n_se_tc2[:,3],markercolor=:red,markersize=5)
scatter!(n_se_tc3[:,1],n_se_tc3[:,2],yerror=n_se_tc3[:,3],markercolor=:cyan,markersize=5,legend=false)
# scatter!(n_se_c[:,1],n_se_c[:,2],yerror=n_se_c[:,3],markercolor=:green,markersize=5,legend=false)
##
plot()
plot!(log10.(fit_t.t),fit_t.n,yerror=fit_t.n_std,c=:black,lw=5)
plot!(log10.(fit_tc1.t),fit_tc1.n,yerror=fit_tc1.n_std,c=:blue,lw=5)
plot!(log10.(fit_tc2.t),fit_tc2.n,yerror=fit_tc2.n_std,c=:red,lw=5)
plot!(log10.(fit_tc3.t),fit_tc3.n,yerror=fit_tc3.n_std,c=:cyan,lw=5)

scatter!(log10.(n_se_t[:,1]),  n_se_t[:,2],  yerror=n_se_t[:,3],  markercolor=:black,markersize=5)
scatter!(log10.(n_se_tc1[:,1]),n_se_tc1[:,2],yerror=n_se_tc1[:,3],markercolor=:blue,markersize=5)
scatter!(log10.(n_se_tc2[:,1]),n_se_tc2[:,2],yerror=n_se_tc2[:,3],markercolor=:red,markersize=5)
scatter!(log10.(n_se_tc3[:,1]),n_se_tc3[:,2],yerror=n_se_tc3[:,3],markercolor=:cyan,markersize=5,legend=false)

##
res_dir = main_dir*"TCR-CD4 modeling/results_adhesion/"

CSV.write(res_dir*"fit_t.csv", fit_t)
CSV.write(res_dir*"fit_tc1.csv", fit_tc1)
CSV.write(res_dir*"fit_tc2.csv", fit_tc2)
CSV.write(res_dir*"fit_tc3.csv", fit_tc3)


##
plot()
plot!(fit_t.t,fit_t.Pa,yerror=fit_t.Pa_std,c=:black,lw=5)
plot!(fit_tc1.t,fit_tc1.Pa,yerror=fit_tc1.Pa_std,c=:blue,lw=5)
plot!(fit_tc2.t,fit_tc2.Pa,yerror=fit_tc2.Pa_std,c=:red,lw=5)
plot!(fit_tc3.t,fit_tc3.Pa,yerror=fit_tc3.Pa_std,c=:cyan,lw=5)


CSV.write(res_dir*"k_tmc_af.csv", rates2Df(k_tc))


# scatter!(n_se_t[:,1],  Pa_t,  markercolor=:black,markersize=5)
# scatter!(n_se_tc1[:,1],Pa_tc1[:,2],markercolor=:blue,markersize=5)
# scatter!(n_se_tc2[:,1],n_se_tc2[:,2],yerror=n_se_tc2[:,3],markercolor=:red,markersize=5)
# scatter!(n_se_tc3[:,1],n_se_tc3[:,2],yerror=n_se_tc3[:,3],markercolor=:cyan,markersize=5,legend=false)

# plot()
# f = (du,u,p,t) -> trimolecular_no_cons!(du,u,p,t,densities_T_full, [])
# Color = :black; fit_t = plot_sum(rates,:black)
# f = (du,u,p,t) -> trimolecular_no_cons!(du,u,p,t,densities_TC1, [])
# Color = :blue; fit_tc1 = plot_sum(rates,:blue)
# # plot_sum([0., 0.])
# f = (du,u,p,t) -> trimolecular_no_cons!(du,u,p,t,densities_TC2, [])
# Color = :red; fit_tc2 = plot_sum(rates,:red)
# # plot_sum([0., 0.])
# f = (du,u,p,t) -> trimolecular_no_cons!(du,u,p,t,densities_TC3, [])
# Color = :cyan; fit_tc3 = plot_sum(rates,:cyan)

# scatter!(n_se_tc1[:,1],n_se_tc1[:,2],yerror=n_se_tc1[:,3],markercolor=:blue,markersize=5)
# scatter!(n_se_tc2[:,1],n_se_tc2[:,2],yerror=n_se_tc2[:,3],markercolor=:red,markersize=5)
# scatter!(n_se_tc3[:,1],n_se_tc3[:,2],yerror=n_se_tc3[:,3],markercolor=:cyan,markersize=5,legend=false)