# TO RUN THIS SCRIPT, FIRST YOU NEED TO HAVE RUN "TCR_CD4_MHC_adhesion_freq_fit.jl" to load all relevant data and dependencies

# from adhesion frequency fitting:
rates_tc = [0.000831216  # AcKaTCR
            0.526597651  # k₋₁
            0.05739801   # KaCD4*
            1.992950078] # k₋₂

densities_TrvMr = [42.0, 5.0, 1.0] # TCR, pMHC, CD4, Ac
densities_TCrvMr = [14.0, 29.0, 11.0, 1.0] # TCR, pMHC, CD4, Ac
densities_TCrvMb = [2.0, 52.0, 2.0, 1.0] # TCR, pMHC, CD4, Ac
densities_TCbvMr = [49.0, 4.0, 53.0, 1.0] # TCR, pMHC, CD4, Ac

TrvMr_data = CSV.read(data_dir*"TCR_RBC_v_MHC_RBC.csv", DataFrame)
TCrvMr_data = CSV.read(data_dir*"TCR_CD4_RBC_v_MHC_RBC.csv", DataFrame)
TCrvMb_data = CSV.read(data_dir*"TCR_CD4_RBC_v_MHC_bead.csv", DataFrame)
TCbvMr_data = CSV.read(data_dir*"TCR_CD4_bead_v_MHC_RBC.csv", DataFrame)


t_TrvMr, n_TrvMr, n_se_TrvMr = clean_n_norm_data(TrvMr_data,  densities_TrvMr)
t_TCrvMr, n_TCrvMr, n_se_TCrvMr = clean_n_norm_data(TCrvMr_data,  densities_TCrvMr)
t_TCrvMb, n_TCrvMb, n_se_TCrvMb = clean_n_norm_data(TCrvMb_data,  densities_TCrvMb)
t_TCbvMr, n_TCbvMr, n_se_TCbvMr = clean_n_norm_data(TCbvMr_data,  densities_TCbvMr)

##
plot()
scatter(n_se_TrvMr[:,1],n_se_TrvMr[:,2],yerror=n_se_TrvMr[:,3],markercolor=:green,markersize=5,label="TCR_RBC_v_MHC_RBC")
scatter!(n_se_TCrvMr[:,1],n_se_TCrvMr[:,2],yerror=n_se_TCrvMr[:,3],markercolor=:black,markersize=5,label="TCR_CD4_RBC_v_MHC_RBC")
scatter!(n_se_TCrvMb[:,1],n_se_TCrvMb[:,2],yerror=n_se_TCrvMb[:,3],markercolor=:blue,markersize=5,label="TCR_CD4_RBC_v_MHC_bead")
scatter!(n_se_TCbvMr[:,1],n_se_TCbvMr[:,2],yerror=n_se_TCbvMr[:,3],markercolor=:orange,markersize=5,label="TCR_CD4_bead_v_MHC_RBC", legend=:outertopright)

#________________________________________________________________________________________________________________________
## FIT TCR_RBC_v_MHC_RBC by bimolecular model

tspan = (0., 16. +1)
u₀ = [0.0]
pinit = rates_tc[1:2] # [k₁, k₋₁]
bounds = [(1.0e-6,100.0), (1.0e-6,100.0)] # [k₁, k₋₁]


TrvMr = Data([n_TrvMr],[t_TrvMr])
# TrvMr = Data([n_t],[t_t])

M_TrvMr = M_global(bimolecular!,u₀,tspan,[densities_TrvMr],[[]],0.0,Rodas5)

res_TrvMr = global_LS_fit(M_TrvMr,TrvMr,bounds)
rates_TrvMr = res_TrvMr.rates
σ_TrvMr = hessian2σ(res_TrvMr.hes,res_TrvMr.l,TrvMr.n[1],rates_TrvMr)
@show k_TrvMr = rates_TrvMr .± σ_TrvMr

fit_TrvMr = plot_n_fit(M_TrvMr,k_TrvMr,densities_TrvMr,[];Color=:green)

rates_TrvMr_conv = [rates_TrvMr[1]/rates_TrvMr[2], rates_TrvMr[2]] # convert to [AcKₐT, k₋₁] from [k₁, k₋₁], with AcKₐT = k₁/k₋₁

# uncomment to save results
# CSV.write(res_dir*"fit_TrvMr.csv", fit_TrvMr)
# CSV.write(res_dir*"k_TrvMr.csv", rates2Df(rates_TrvMr_conv))

#________________________________________________________________________________________________________________________
## fit TCR_RBC_v_MHC_RBC and "TCR_CD4_RBC_v_MHC_RBC by trimolecular model and global fitting

## KEEP f_tol to 1e-5
u₀ = [0.0, 0.0] # [n_TCR, n_3]
cons_TCrvMr = fill([],4) # don't hold any rate constant
densities_TrvMr_full = [42.0, 5.0, 0.0, 1.0]
dens_TCrvMr = [densities_TrvMr_full, densities_TCrvMr]
tspan = (0,16+1)
# pinit = [0.0,0.0,0.0,0.0] # [k₁, k₋₁, k₂, k₋₂]
pinit = rates_tc


TCrvMr = Data([n_TrvMr, n_TCrvMr],[t_TrvMr, t_TCrvMr]) # fit all experiments at same time
model = trimolecular_ana_zhu_AcKa

M_TCrvMr = M_global(model,u₀,tspan,dens_TCrvMr,cons_TCrvMr,0.0,Rodas5)

bounds = [(1.0e-7,10.0),(1.0e-7,10.0), # [AcKaTCR, k₋₁]
        (1.0e-7, 10),(1.0e-5,100)] # [KaCD4, k₋₂]

res_TCrvMr = global_LS_fit(M_TCrvMr,TCrvMr,bounds)
rates_TCrvMr = res_TCrvMr.rates
res_TCrvMr.hes
σ_TCrvMr = hessian2σ(res_TCrvMr.hes,res_TCrvMr.l,reduce(vcat,TCrvMr.n),res_TCrvMr.rates) #./ √(64)
k_TCrvMr = rates_TCrvMr .± σ_TCrvMr

plot_n_fit(M_TCrvMr,k_TCrvMr,densities_TrvMr_full,[];Color=:green)
fit_TCrvMr = plot_n_fit(M_TCrvMr,k_TCrvMr,densities_TCrvMr,    [];Color=:black)

# uncomment to save results
# CSV.write(res_dir*"fit_TCrvMr.csv", fit_TCrvMr)
# CSV.write(res_dir*"k_TCrvMr.csv", rates2Df(k_TCrvMr))



#________________________________________________________________________________________________________________________
# try fitting the TCR_CD4_RBC_v_MHC_RBC interaction by keeping the bimolecular TCR rates constant
# DOES NOT FIT WELL, (cuts through it, visually)

# model = trimolecular_ana_zhu_cons_bi
# u₀ = [0.0, 0.0] # [n_TCR, n_3]
# dens_TCrvMr = [densities_TCrvMr]
# cons_TCrvMr = rates_TrvMr_conv # don't hold any rate constant
# TCrvMr = Data([n_TCrvMr],[t_TCrvMr]) # fit all experiments at same time
# M_TCrvMr = M_global(model,u₀,tspan,dens_TCrvMr,[cons_TCrvMr],0.0,Rodas5)
# bounds = [(1.0e-8, 10),(1.0e-8,10)] # [KaCD4, k₋₂]
# pinit = rates_tc[3:4]
# res_TCrvMr = global_LS_fit(M_TCrvMr,TCrvMr,bounds)
# rates_TCrvMr = res_TCrvMr.rates
# res_TCrvMr.hes
# σ_TCrvMr = hessian2σ(res_TCrvMr.hes,res_TCrvMr.l,reduce(vcat,TCrvMr.n),res_TCrvMr.rates) #./ √(64)
# k_TCrvMr = rates_TCrvMr .± σ_TCrvMr

# plot_n_fit(M_TCrvMr,k_TCrvMr,densities_TCrvMr,cons_TCrvMr;Color=:black)
# res_TCrvMr.l


#________________________________________________________________________________________________________________________
## fit all trimolecular interactions by trimolecular model
model = trimolecular_ana_zhu_AcKa
u₀ = [0.0, 0.0] # [n_TCR, n_3]
tspan = (0,16+1)
pinit = rates_tc

#________________________________________________________________________________________________________________________
## fit TCR_CD4_RBC_v_MHC_RBC by trimolecular model
# DOES not fit well (very large uncertainty, it's overfitting, as we are fitting 4 params with little data, 14 pts)

#= #uncomment here
dens_TCrvMr = [densities_TCrvMr]
TCrvMr = Data([n_TCrvMr],[t_TCrvMr]) # fit all experiments at same time
M_TCrvMr = M_global(model,u₀,tspan,dens_TCrvMr,[[]],0.0,Rodas5)
bounds = [(1.0e-6,10.0),(1.0e-9,10.0), # [AcKaTCR, k₋₁]
        (1.0e-6, 10),(1.0e-5,10)] # [KaCD4, k₋₂]
# bounds = [(0,Inf),(0,Inf), # [AcKaTCR, k₋₁]
# (       0,Inf),(0,Inf)] # [KaCD4, k₋₂]
res_TCrvMr = global_LS_fit(M_TCrvMr,TCrvMr,bounds)
rates_TCrvMr = res_TCrvMr.rates
res_TCrvMr.hes
σ_TCrvMr = hessian2σ(res_TCrvMr.hes,res_TCrvMr.l,reduce(vcat,TCrvMr.n),res_TCrvMr.rates) #./ √(64)
k_TCrvMr = rates_TCrvMr .± σ_TCrvMr

plot_n_fit(M_TCrvMr,k_TCrvMr,densities_TCrvMr,[];Color=:black)
=#

#________________________________________________________________________________________________________________________
## fit TCR_CD4_RBC_v_MHC_bead
# DOES not fit well (very large uncertainty, it's overfitting, as we are fitting 4 params with little data)

#=
dens_TCrvMb = [densities_TCrvMb]
TCrvMb = Data([n_TCrvMb],[t_TCrvMb]) # fit all experiments at same time
M_TCrvMb = M_global(model,u₀,tspan,dens_TCrvMb,[[]],0.0,Rodas5)
bounds = [(1.0e-6,10.0),(1.0e-9,10.0), # [AcKaTCR, k₋₁]
        (1.0e-6, 10),(1.0e-5,10)] # [KaCD4, k₋₂]
pinit = rates_tc

res_TCrvMb = global_LS_fit(M_TCrvMb,TCrvMb,bounds)
rates_TCrvMb = res_TCrvMb.rates
res_TCrvMb.hes
σ_TCrvMb = hessian2σ(res_TCrvMb.hes,res_TCrvMb.l,reduce(vcat,TCrvMb.n),res_TCrvMb.rates) #./ √(64)
k_TCrvMb = rates_TCrvMb .± σ_TCrvMb

scatter(n_se_TCrvMb[:,1],n_se_TCrvMb[:,2],yerror=n_se_TCrvMb[:,3],markercolor=:blue,markersize=5,label="TCR_CD4_RBC_v_MHC_bead")
plot_n_fit(M_TCrvMb,k_TCrvMb,densities_TCrvMb,[];Color=:blue)
=#

#________________________________________________________________________________________________________________________
## fit TCR_CD4_bead_v_MHC_RBC by trimolecular model
# 
## KEEP f_tol to 1e-5 or global_LS_fit in model_fit_global_v2.jl file

# #=

dens_TCbvMr = [densities_TCbvMr]
TCbvMr = Data([n_TCbvMr],[t_TCbvMr]) # fit all experiments at same time
M_TCbvMr = M_global(model,u₀,tspan,dens_TCbvMr,[[]],0.0,Rodas5)
bounds = [(1.0e-6,10.0),(1.0e-9,10.0), # [AcKaTCR, k₋₁]
        (1.0e-9, 10),(1.0e-5,10)] # [KaCD4, k₋₂]
res_TCbvMr = global_LS_fit(M_TCbvMr,TCbvMr,bounds)
rates_TCbvMr = res_TCbvMr.rates
res_TCbvMr.hes
σ_TCbvMr = hessian2σ(res_TCbvMr.hes,res_TCbvMr.l,reduce(vcat,TCbvMr.n),res_TCbvMr.rates) #./ √(64)
k_TCbvMr = rates_TCbvMr .± σ_TCbvMr

fit_TCbvMr = plot_n_fit(M_TCbvMr,k_TCbvMr,densities_TCbvMr,[];Color=:orange)
scatter(n_se_TCbvMr[:,1],n_se_TCbvMr[:,2],yerror=n_se_TCbvMr[:,3],markercolor=:orange,markersize=5,label="TCR_CD4_bead_v_MHC_RBC", legend=:outertopright)

# uncomment to save
# CSV.write(res_dir*"fit_TCbvMr.csv", fit_TCbvMr)
# CSV.write(res_dir*"k_TCbvMr.csv", rates2Df(k_TCbvMr))

# =#

## since the two bead reactions are way overfitting (really large standard errors), try fitting them by keeping TCR on and off rates fixed.
model = trimolecular_ana_zhu_cons_k₋₁_k₋₂
#________________________________________________________________________________________________________________________
## fit TCR_CD4_RBC_v_MHC_bead, keeping off rates constant
# THIS FITS WELL
dens_TCrvMb = [densities_TCrvMb]
cons_from_tc = rates_tc[[2,4]] # keep [k₋₁, k₋₂] costant
cons_from_tc = rates_TCbvMr[[2,4]] 

TCrvMb = Data([n_TCrvMb],[t_TCrvMb]) # fit all experiments at same time
M_TCrvMb = M_global(model,u₀,tspan,dens_TCrvMb,[cons_from_tc],0.0,Rodas5)
bounds = [(1.0e-6,10.0),
            (1.0e-6, 10)] # [KaCD4, k₋₂]
pinit = rates_tc[[1,3]]
res_TCrvMb = global_LS_fit(M_TCrvMb,TCrvMb,bounds)
rates_TCrvMb = res_TCrvMb.rates
res_TCrvMb.hes
σ_TCrvMb = hessian2σ(res_TCrvMb.hes,res_TCrvMb.l,reduce(vcat,TCrvMb.n),res_TCrvMb.rates) #./ √(64)
k_TCrvMb = rates_TCrvMb .± σ_TCrvMb

plot_n_fit(M_TCrvMb,k_TCrvMb,densities_TCrvMb,cons_from_tc;Color=:blue)

#________________________________________________________________________________________________________________________
## fit TCR_CD4_RBC_v_MHC_bead with bi rates constant

model = trimolecular_ana_zhu_cons_bi
dens_TCrvMb = [densities_TCrvMb]
# cons_from_tc = rates_tc[[1,2]] # keep [k₋₁, k₋₂] costant
cons_from_tc = rates_TCbvMr[[1,2]]

TCrvMb = Data([n_TCrvMb],[t_TCrvMb]) # fit all experiments at same time
M_TCrvMb = M_global(model,u₀,tspan,dens_TCrvMb,[cons_from_tc],0.0,Rodas5)
bounds = [(1.0e-9,10.0),
            (1.0e-9, 10)] # [KaCD4, k₋₂]
# pinit = rates_tc[[3,4]]
pinit = rates_TCbvMr[[3,4]]
res_TCrvMb = global_LS_fit(M_TCrvMb,TCrvMb,bounds)
rates_TCrvMb = res_TCrvMb.rates
res_TCrvMb.hes
σ_TCrvMb = hessian2σ(res_TCrvMb.hes,res_TCrvMb.l,reduce(vcat,TCrvMb.n),res_TCrvMb.rates) #./ √(64)
k_TCrvMb = rates_TCrvMb .± σ_TCrvMb

plot_n_fit(M_TCrvMb,k_TCrvMb,densities_TCrvMb,cons_from_tc;Color=:blue)

#________________________________________________________________________________________________________________________
## fit TCR_CD4_RBC_v_MHC_bead with k₋₁ rate constant (USE THIS FOR PAPER FIG 2H)

# KEEP f_tol low at 1e-8 !!!
model = trimolecular_ana_zhu_cons_k₋₁

cons_k₋₁ = rates_TCbvMr[2]
# cons_k₋₁ = rates_tc[2]
# cons_k₋₁ = rates_TCrvMr[2]

M_TCrvMb = M_global(model,u₀,tspan,dens_TCrvMb,cons_k₋₁,0.0,Rodas5)
bounds = [(1.0e-8,10.0),
          (0, 10), (1.0e-8,10.0)]
pinit = rates_tc[[1,3,4]] #  [KaTCR, KaCD4, k₋₂]


res_TCrvMb = global_LS_fit(M_TCrvMb,TCrvMb,bounds)
rates_TCrvMb = res_TCrvMb.rates
res_TCrvMb.hes
σ_TCrvMb = hessian2σ(res_TCrvMb.hes,res_TCrvMb.l,reduce(vcat,TCrvMb.n),res_TCrvMb.rates) #./ √(64)
k_TCrvMb = rates_TCrvMb .± σ_TCrvMb #  [KaTCR, KaCD4, k₋₂]

fit_TCrvMb = plot_n_fit(M_TCrvMb,k_TCrvMb,densities_TCrvMb,cons_k₋₁;Color=:blue)

# uncomment for saving
# CSV.write(res_dir*"fit_TCrvMb.csv", fit_TCrvMb)
# CSV.write(res_dir*"k_TCrvMb.csv", rates2Df(k_TCrvMb))


#________________________________________________________________________________________________________________________
# FITTING OF OTHER ALTERNATIVE MODELS CONSIDERED, NOT INCLUDED IN FIGURE 

## fit TCR_CD4_bead v MHC_RBC by trimolecular model (orange)
# model = trimolecular_ana_zhu_cons_k₋₁_k₋₂

# cons_from_tc = rates_tc[[2,4]] # keep [k₋₁, k₋₂] costant
# cons_from_tc = rates_TCrvMr[[2,4]]
# dens_TCbvMr = [densities_TCbvMr]
# TCbvMr = Data([n_TCbvMr],[t_TCbvMr]) # fit all experiments at same time
# M_TCbvMr = M_global(model,u₀,tspan,dens_TCbvMr,[cons_from_tc],0.0,Rodas5)
# bounds = [(1.0e-6,1.0),
#         (0, 1)]
# pinit = rates_tc[[1,3]]
# res_TCbvMr = global_LS_fit(M_TCbvMr,TCbvMr,bounds)
# rates_TCbvMr = res_TCbvMr.rates
# res_TCbvMr.hes
# σ_TCbvMr = hessian2σ(res_TCbvMr.hes,res_TCbvMr.l,reduce(vcat,TCbvMr.n),res_TCbvMr.rates) #./ √(64)
# # σ_TCbvMr[2] = 0
# k_TCbvMr = rates_TCbvMr .± σ_TCbvMr


# plot_n_fit(M_TCbvMr,k_TCbvMr,densities_TCbvMr,cons_from_tc;Color=:orange)
# scatter(n_se_TCbvMr[:,1],n_se_TCbvMr[:,2],yerror=n_se_TCbvMr[:,3],markercolor=:orange,markersize=5,label="TCR_CD4_bead_v_MHC_RBC", legend=:outertopright)


# model = trimolecular_ana_zhu_cons_k₋₁
# cons_k₋₁ = rates_tc[2]
# cons_k₋₁ = 10.7
# M_TCbvMr = M_global(model,u₀,tspan,dens_TCbvMr,cons_k₋₁,0.0,Rodas5)
# bounds = [(1.0e-8,10.0),
#         (0, 10), (1.0e-8,10.0)]
# pinit = rates_tc[[1,3,4]] #  [KaTCR, KaCD4, k₋₂]
# pinit[2] = 0

# res_TCbvMr = global_LS_fit(M_TCbvMr,TCbvMr,bounds)
# rates_TCbvMr = res_TCbvMr.rates
# res_TCbvMr.hes
# σ_TCbvMr = hessian2σ(res_TCbvMr.hes,res_TCbvMr.l,reduce(vcat,TCbvMr.n),res_TCbvMr.rates) #./ √(64)
# k_TCbvMr = rates_TCbvMr .± σ_TCbvMr #  [KaTCR, KaCD4, k₋₂]

# plot_n_fit(M_TCbvMr,k_TCbvMr,densities_TCbvMr,cons_k₋₁;Color=:orange)


# model = trimolecular_ana_zhu_cons_bi
# dens_TCbvMr = [densities_TCbvMr]
# cons_from_tc = rates_tc[[1,2]] # keep [k₋₁, k₋₂] costant
# cons_from_tc = rates_TCrvMr[[1,2]]

# M_TCbvMr = M_global(model,u₀,tspan,dens_TCbvMr,[cons_from_tc],0.0,Rodas5)
# bounds = [(0, 10), (1.0e-8,100.0)]
# pinit = [1e-8,10] #  [KaTCR, KaCD4, k₋₂]

# res_TCbvMr = global_LS_fit(M_TCbvMr,TCbvMr,bounds)
# rates_TCbvMr = res_TCbvMr.rates
# res_TCbvMr.hes
# σ_TCbvMr = hessian2σ(res_TCbvMr.hes,res_TCbvMr.l,reduce(vcat,TCbvMr.n),res_TCbvMr.rates) #./ √(64)
# k_TCbvMr = rates_TCbvMr .± σ_TCbvMr #  [KaTCR, KaCD4, k₋₂]

# plot_n_fit(M_TCbvMr,rates_TCbvMr,densities_TCbvMr,cons_from_tc;Color=:orang