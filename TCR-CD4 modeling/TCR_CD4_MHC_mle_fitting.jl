using Revise
using StatsBase, KernelDensity, Distributions, Measurements
using DataFrames, Plots, ColorSchemes
main_dir = "C:/Users/stravaglino3/Documents/TCR-MHC-CD4 modeling for github/"

src_dir = main_dir*"src/"
includet(src_dir*"models.jl")
includet(src_dir*"models_dissociation.jl")
includet(src_dir*"analytical_models.jl")
includet(src_dir*"helper_functions.jl")
includet(src_dir*"binning_helper_fs.jl")
includet(src_dir*"mle_fitting.jl")
includet(src_dir*"mle_plotting.jl")
using BlackBoxOptim

res_dir = main_dir*"TCR-CD4 modeling/results/"

#______________________________________________________________
##
# load data
data_dir = main_dir*"TCR-CD4 modeling/dissociation data/"

# load bead data
# new bead-bead data from Kaitao combined with Muaz data
tm_b = tm_b_no_boot = CSV.read(data_dir*"bead_TCR_MHC_combined_no_outl_new.csv", DataFrame) #combined KL MR
cm_b = cm_b_no_boot = CSV.read(data_dir*"CD4_LT_bead_KL.csv", DataFrame) # here we use Kaitao's data only
tmc_b = tmc_b_no_boot = CSV.read(data_dir*"bead_TCR_MHC_CD4_combined.csv", DataFrame); #combined KL MR

thermal_b = CSV.read(data_dir*"Thermal fluctuation assay _ CD4-E8TCR-TPI_HLADR.csv",DataFrame)
tm_b_thrml = delete_missing(thermal_b[:,3:4])[1]
cm_b_thrml = delete_missing(thermal_b[:,8:9])[1]
tmc_b_thrml = delete_missing(thermal_b[:,13:14])[1]

# tm_b = combined_bootstrap_dataset([tm_b.F,tm_b.t], 30)
# # cm_b = combined_bootstrap_dataset(cm_b, 20)
# tmc_b = combined_bootstrap_dataset([tmc_b.F,tmc_b.t], 30)

# combine 0 force thermal fluctuation data with >0 F-dependent dataset
tm_b_combined = [vcat(tm_b.F,zeros(size(tm_b_thrml))), vcat(tm_b.t,tm_b_thrml)]
cm_b_combined = [vcat(cm_b.F,zeros(size(cm_b_thrml))), vcat(cm_b.t,cm_b_thrml)]
tmc_b_combined = [vcat(tmc_b.F,zeros(size(tmc_b_thrml))), vcat(tmc_b.t,tmc_b_thrml)] 

i_outl = findmax(tmc_b_combined[2])[2] # delete very high outlier
tmc_b_combined[1] = deleteat!(tmc_b_combined[1], i_outl)
tmc_b_combined[2] = deleteat!(tmc_b_combined[2], i_outl)



tm_b_combined = combined_bootstrap_dataset(tm_b_combined, 30)
# cm_b = combined_bootstrap_dataset(cm_b, 20)
tmc_b_combined = combined_bootstrap_dataset(tmc_b_combined, 30)

# tm_b = F_t_Data(tm_b.F,tm_b.t,   InterpKDE(kde(tm_b.F)))
tm_b = F_t_Data(tm_b_combined[1],tm_b_combined[2],  InterpKDE(kde(tm_b_combined[1])))
cm_b  = F_t_Data(cm_b_combined[1],cm_b_combined[2],   InterpKDE(kde(cm_b_combined[1])))
tmc_b = F_t_Data(tmc_b_combined[1],tmc_b_combined[2], InterpKDE(kde(tmc_b_combined[1])))

#______________________________________________________________
## Plotting dataset
gr()
s1=Plots.scatter(tm_b.F, tm_b.t, label="TCR-MHC");
s2=Plots.scatter(cm_b.F, cm_b.t, label="MHC-CD4");
s3=Plots.scatter(tmc_b.F, tmc_b.t, label="TCR-MHC-CD4");
plot(s1,s2,s3)

h1=histogram2d(tm_b.F, log10.(tm_b.t), nbins=(40,50),title="TCR-MHC")
h2=histogram2d(cm_b.F, log10.(cm_b.t), nbins=(40,50),title="MHC-CD4")
h3=histogram2d(tmc_b.F, log10.(tmc_b.t), nbins=(40,50),title="TCR-MHC-CD4")#, c=:turbo);#,aspect_ratio=2)
plot(h1,h2,h3)

#______________________________________________________________
## BINNING

# edges=[0,11,19,23,30,33.,45]
edges=[0,11,19,22.5,29,33.,45]
# edges=nothing
TM_B = bootstrap_bins([tm_b_no_boot.F, tm_b_no_boot.t],6,300;bin_type="force",edges=edges) # [ΓΕ.F, ΓΕ.t]
plot()
clr=:blue; p1 = plot_F_vs_t(TM_B)
CM_B = bootstrap_bins([cm_b_no_boot.F, cm_b_no_boot.t],6,300;bin_type="force",edges=nothing) # [ΓΕ.F, ΓΕ.t]
clr=:red; p1 = plot_F_vs_t(CM_B)
edges=[0,8,13,17.5,22.5,27,33,45]
TMC_B = bootstrap_bins([tmc_b_no_boot.F, tmc_b_no_boot.t],8,300;bin_type="force",edges=edges) # [ΓΕ.F, ΓΕ.t]
clr=:green; p1 = plot_F_vs_t(TMC_B)
# plot!(0:40,mean_t_tmc_b,lw=5,c=:green,label="TCR-MHC-CD4")


#_________________________________________________________
## fit TCR-MHC catch bond (letting kf free to change)
u₀ = [] # don't provide initial states
u₀_type = "equilibrium" # calculate initial states based on equilibrium
# u₀_type = "detailed_balance" # calculate initial states based on equilibrium
tspan = (0.0,50)
# kf_tm = 1/mean(tm_b_thrml)
model = bimolecular_diss_bell!
M_tm_b  = M(model,u₀,u₀_type,tspan,[],[],0.0,Rodas5)

pinit = [ 10, 0.01,
          0.1, 0.1,
          1.0, 0.1,
        100.0, 1.0]

bounds = [ (10,100), (0.,10.),
        (0.,10.),(0.,10.),
        (0.,10.),(0.,10.),
        (0.,100.),(0.,10.)]

@time mle_loss(pinit, tm_b, M_tm_b)

model = two_state_catch_bond_analytical
M_tm_b  = M(model,u₀,u₀_type,tspan,[],[],0.0,Rodas5)
@time mle_loss(pinit, tm_b, M_tm_b)

using Random
Random.seed!(3)

@time res_tm_b = mle_fit(M_tm_b,tm_b,bounds)
ks_tm_b = res_tm_b.rates
σ_tm_b = hessian2σ(res_tm_b.hes)
@show k_tm_b = ks_tm_b .± σ_tm_b

u₀tm_b = initial_u(ks_tm_b,M_tm_b)
# CSV.write(res_dir*"k_γϵ.csv", rates2Df(k_γϵ))

# plotly()
plot()
tm_b_M_w_uncert = M(model,u₀,u₀_type,tspan,[],[],0.0,Rodas5)
mean_t_tm_b = mean_ft_fit_analytical_2(tm_b_M_w_uncert,k_tm_b,tm_b,F_max=40.,t_max=100.)
mean_t_tm_b = map(t->Measurement(t), mean_t_tm_b)
plot!(0:40,mean_t_tm_b,lw=5,c=:blue,label="TCR-MHC")
clr = :blue
plot_F_vs_t(TM_B)

# tm_b_M_w_uncert = M(two_state_catch_bond_ana_matlab,u₀,u₀_type,tspan,[],[],0.0,Rodas5)
# mean_t_tm_b = mean_ft_fit_analytical_2(tm_b_M_w_uncert,ks_tm_b,tm_b,F_max=40.,t_max=100.)
# mean_t_tm_b = map(t->Measurement(t), mean_t_tm_b)
# plot!(0:40,mean_t_tm_b,lw=5,c=:red,linestyle=:dash,label="TCR-MHC")

tm_b_M_w_uncert = M(bimolecular_diss_bell!,u₀,u₀_type,tspan,[],[],0.0,Rodas5)
plot_mean_ft_fit(tm_b_M_w_uncert,ks_tm_b,tm_b,TM_B;species="TCR-MHC")

#_____________________________________________________
## fit CD4-MHC slip bond

u₀ = [1.0] # don't provide initial states
u₀_type = [] # calculate initial states based on equilibrium
# u₀_type = "detailed_balance" # calculate initial states based on equilibrium
tspan = (0.0,50)

model = bimolecular_diss_bell_single!
M_cm_b  = M(model,u₀,u₀_type,tspan,[],[],0.0,Rodas5)

pinit = [1.0, 0.01]

bounds = [(0.,10.),(0.,10.)]

@time mle_loss(pinit, cm_b, M_cm_b)

model = slip_bond_analytical
M_cm_b  = M(model,u₀,u₀_type,tspan,[],[],0.0,Rodas5)
@time mle_loss(pinit, cm_b, M_cm_b)

@time res_cm_b = mle_fit(M_cm_b,cm_b,bounds)
ks_cm_b = res_cm_b.rates
σ_cm_b = hessian2σ(res_cm_b.hes)
@show k_cm_b = ks_cm_b .± σ_cm_b

u₀cm_b = initial_u(ks_cm_b,M_cm_b)

plot()
cm_b_M_w_uncert = M(model,u₀,u₀_type,tspan,[],[],0.0,Rodas5)
mean_t_cm_b = mean_ft_fit_analytical_2(cm_b_M_w_uncert,k_cm_b,cm_b,F_max=40.,t_max=100.)
mean_t_cm_b = map(t->Measurement(t), mean_t_cm_b)
plot!(0:40,mean_t_cm_b,lw=5,c=:red,label="CD4-MHC")
clr = :red
plot_F_vs_t(CM_B)



#_____________________________________________________
## trimolecular fit

u₀ = [] # don't provide initial states
u₀_type = "equilibrium" # calculate initial states based on equilibrium
# u₀_type = "detailed_balance" # calculate initial states based on equilibrium
tspan = (0.0,20)
model = two_state_catch_bond_analytical
M_tmc_b  = M(model,u₀,u₀_type,tspan,[],[],0.0,Rodas5)


pinit = [ 6, 0.,
        0.05,0.4,
        4,0.05,
        30,0.8]

bounds = [(1.,10.),(0.,10),
        (0.,0.05),(0.,10.),
        (0.,10.),(0.,10),
        (0.,100.),(0.,10)]
@time mle_loss(pinit, tmc_b, M_tmc_b)

@time res_tmc_b = mle_fit(M_tmc_b,tmc_b,bounds)
ks_tmc_b = res_tmc_b.rates
σ_tmc_b = hessian2σ(res_tmc_b.hes)
@show k_tmc_b = ks_tmc_b .± σ_tmc_b

        
# ks_tmc_b[3] = 0.045
tmc_b_M_w_uncert = M(model,u₀,u₀_type,tspan,[],ks_tm_b,0.0,Rodas5)
mean_t_tmc_b = mean_ft_fit_analytical_2(tmc_b_M_w_uncert,k_tmc_b,tmc_b,F_max=40.,t_max=100.)
# mean_t_tmc_b = mean_ft_fit_analytical_2(tmc_b_M_w_uncert,pinit,tmc_b,F_max=40.,t_max=100.)

plot()
mean_t_tmc_b = map(t->Measurement(t), mean_t_tmc_b)
plot!(0:40,mean_t_tmc_b,lw=5,c=:green,label="TCR-MHC-CD4")
clr = :green
plot_F_vs_t(TMC_B)


tmc_b_M_num = M(bimolecular_diss_bell!,u₀,u₀_type,tspan,[],ks_tm_b,0.0,Rodas5)
plot_mean_ft_fit(tmc_b_M_num,k_tmc_b,tmc_b,TMC_B;species="TCR-MHC-CD4")

# plot()
# tmc_b_M_w_uncert = M(model,u₀,u₀_type,tspan,[],ks_tm_b,0.0,Rodas5)
# mean_t_tmc_b = mean_ft_fit_analytical_2(tmc_b_M_w_uncert,k_tmc_b,tmc_b,F_max=40.,t_max=100.)
# mean_t_tmc_b = map(t->Measurement(t), mean_t_tmc_b)
# plot!(0:40,mean_t_tmc_b,lw=5,c=:blue,label="TCR-MHC")
# clr = :blue
# plot_F_vs_t(TMC_B)

# # tm_b_M_w_uncert = M(two_state_catch_bond_ana_matlab,u₀,u₀_type,tspan,[],[],0.0,Rodas5)
# # mean_t_tm_b = mean_ft_fit_analytical_2(tm_b_M_w_uncert,ks_tm_b,tm_b,F_max=40.,t_max=100.)
# # mean_t_tm_b = map(t->Measurement(t), mean_t_tm_b)
# # plot!(0:40,mean_t_tm_b,lw=5,c=:red,linestyle=:dash,label="TCR-MHC")

# tmc_b_M_num = M(bimolecular_diss_bell!,u₀,u₀_type,tspan,[],ks_tm_b,0.0,Rodas5)
# plot_mean_ft_fit(tmc_b_M_num,ks_tmc_b,tmc_b,TMC_B;species="TCR-MHC")

##______________________________________________________________
# plotting for paper

plot()
rates_tm_b_by_F = force_dependence_of_catch_rates(ks_tm_b); ylims!((0,20)); title!("TCR-MHC");
rates_tmc_b_by_F = force_dependence_of_catch_rates(ks_tmc_b); ylims!((0,20)); title!("TCR-MHC-CD4");
rates_cm_b_by_F = force_dependence_of_catch_rates(ks_cm_b); ylims!((0,50)); title!("CD4-MHC");


function states_proportion(ka,k₋a)
        (10 ^ ka / (10 ^ ka + 10 ^k₋a))
end

s_w_tm_b = states_proportion.(rates_tm_b_by_F.ka, rates_tm_b_by_F.k₋a)
s_w_tmc_b = states_proportion.(rates_tmc_b_by_F.ka, rates_tmc_b_by_F.k₋a)
plot(0:0.1:40,s_w_tm_b,lw=5, label="ka/(ka + k₋a) (TCR-MHC)")
plot!(0:0.1:40,s_w_tmc_b,lw=5, label="ka/(ka + k₋a) (TCR-MHC-CD4)")
ylabel!("log₁₀(ka/k₋a)")
xlabel!("F (pN)")

## ----------------------------------------------------------------------------



# LIFETIME COOPERATIVITY FUNCTIONS
function plot_lftm_coop(t_bi, t_mix)
        t_bi_pred = t_bi
        # p_1 = plot!(0:35,t_bi_pred,c=:orange,lw=5)
        # a = 1+1
        lftm_coop = lifetime_cooperativity.(t_bi_pred,t_mix)
        plt = plot()
        plt = plot(0:length(t_bi)-1, lftm_coop, errorbar=:ribbon, xlabel="F",ylabel="lifetime cooperativity (%)")
        return lftm_coop, plt
end

function lifetime_cooperativity(t_bi_pred,t_mix)
        lftm_coop = (t_mix - t_bi_pred) / t_bi_pred
end

##-----------------------
# LIFETIME COOPERATIVITY PLOTS
lftm_coop, plt1 = plot_lftm_coop(mean_t_tm_b ,mean_t_tmc_b)
plot(plt1)

plot()
# edges = 0:6:45
edges = [0,6,12,18,24,30,34.5,45]
# edges= nothing
clr = :blue
t_tm_b_same_bins = bootstrap_bins([tm_b_no_boot.F,tm_b_no_boot.t],6,300;bin_type="force",edges=edges) # [ΓΕ.F, ΓΕ.t]
plot_F_vs_t(t_tm_b_same_bins)
clr = :green
t_tmc_b_same_bins = bootstrap_bins([tmc_b_no_boot.F,tmc_b_no_boot.t], 6, 300; bin_type ="force", edges=edges) # [Mix.F, Mix.t]
plot_F_vs_t(t_tmc_b_same_bins)

# lifetime cooperativity at 0 force
t_tm_b_thrml = mean(tm_b_thrml) ± sem(tm_b_thrml);
t_tmc_b_thrml = mean(tmc_b_thrml) ± sem(tmc_b_thrml);
thrml_coop = lifetime_cooperativity(t_tm_b_thrml, t_tmc_b_thrml)
lftm_coop_bins = lifetime_cooperativity.(t_tm_b_same_bins.t_mean .± t_tm_b_same_bins.t_sem,t_tmc_b_same_bins.t_mean .± t_tmc_b_same_bins.t_sem)
plot(plt1)
scatter!(t_tm_b_same_bins.f_mean,lftm_coop_bins)
scatter!([0.],[thrml_coop])
lftm_coop_data = measurements2Df([0.;t_tm_b_same_bins.f_mean],[thrml_coop; lftm_coop_bins]; labels=["lftm coop (%)","std"])

CSV.write(res_dir*"t_tm_b_fit.csv", measurements2Df(0:40,mean_t_tm_b;labels=["fit","std"]))
CSV.write(res_dir*"t_cm_b_fit.csv", measurements2Df(0:40,mean_t_cm_b;labels=["fit","std"]))
CSV.write(res_dir*"t_tmc_b_fit.csv", measurements2Df(0:40,mean_t_tmc_b;labels=["fit","std"]))
CSV.write(res_dir*"lftm_coop_b.csv", measurements2Df(0:40,lftm_coop;labels=["lftm coop (%)","std"]))

CSV.write(res_dir*"ft_mean_tm_b.csv", TM_B)
CSV.write(res_dir*"ft_mean_cm_b.csv", CM_B)
CSV.write(res_dir*"ft_mean_tmc_b.csv", TMC_B)
CSV.write(res_dir*"lftm_coop_data_b.csv", lftm_coop_data)

CSV.write(res_dir*"bell_rates_tm_b_.csv",rates_tm_b_by_F) 
CSV.write(res_dir*"bell_rates_cm_b_.csv",rates_cm_b_by_F) 
CSV.write(res_dir*"bell_rates_tmc_b_.csv",rates_tmc_b_by_F) 

CSV.write(res_dir*"k_tm_b.csv", rates2Df(k_tm_b))
CSV.write(res_dir*"k_cm_b.csv", rates2Df(k_cm_b))
CSV.write(res_dir*"k_tmc_b.csv", rates2Df(k_tmc_b))



figures_dir = main_dir*"TCR-CD4 modeling/figures/"
pyplot()
plot()
# Plots.scalefontsizes(1/2)
# Plots.scalefontsizes(2)
p1 = plot_survival_surface_analytical(M_tm_b,ks_tm_b,[tm_b_no_boot.F,tm_b_no_boot.t];F_max=45,t_max=10,dF=0.5,dt=0.01)
# savefig(figures_dir*"survival_surface_tm_b.png");  savefig(figures_dir*"survival_surface_tm_b.svg")

plot()
p2 = plot_survival_surface_analytical(M_cm_b,ks_cm_b,[cm_b_no_boot.F,cm_b_no_boot.t];F_max=45,t_max=5,dF=0.5,dt=0.01)
# savefig(figures_dir*"survival_surface_cm_b.png"); savefig(figures_dir*"survival_surface_cm_b.svg")

plot()
p3 = plot_survival_surface_analytical(M_tmc_b,ks_tmc_b,[tmc_b_no_boot.F,tmc_b_no_boot.t];F_max=45,t_max=10,dF=0.5,dt=0.01)
# savefig(figures_dir*"survival_surface_tmc_b.png"); savefig(figures_dir*"survival_surface_tmc_b.svg")



function plot_mean_of_each_species(M,k,data,binned;species=[],dt=0.1,dF=1.0)

        ts = 0:dt:M.tspan[2]
        Fs = 0:dF:maximum(data.F)
        # Fs = 20:10:40

        # Fs = 0:1
    
        M.solver = Tsit5  #change solver to Tsit5 (can handle ± std)
        sols = ftsolve(M,k,Fs,ts)
        _mean_t = mean_lifetime_each_species.(sols)
        # mean_t, proportion = mean_lifetime_each_species.(sols)
        # @show mean_t
        
        mean_t = reduce(hcat,_mean_t)'
        # proportion_stacked = reduce(hcat,proportion)'
        P_activated = mean_t[:,2] ./ (mean_t[:,1] .+ mean_t[:,2])

        plt = plot!(Fs,mean_t,c=clr,
                        # ribbon=0.1,
                        markerstrokecolor=:auto, lw=5,
                        linealpha = 1,
                        ylabel="<t>",
                        xlabel="F",
                        label=species)

                        

        plt = plot!(plot!(Fs,P_activated,c=:black,
        # ribbon=0.1,
                        markerstrokecolor=:auto, lw=5,
                        linealpha = 1,
                        ylabel="<t>",
                        xlabel="F",
                        label="proportion activated",legend=:outerright))
        
        return mean_t, P_activated, plt
end
    
function mean_lifetime_each_species(sol)
        t = sol.t
        # u = sum(sol,dims=1)'
        _u = Array(sol)
        # @show _u[:,1]
        _du = sol(t, Val{1}) 
        # @show size(du)
        # _u = _u ./ _u[:,1]
        
        u = [_u[i,:] for i in 1:size(_u)[1]]
        du = [_du[i,:] for i in 1:size(_du)[1]]

        tᵢ_mean = map(uᵢ -> sum(integrate(t,uᵢ)), u)
        tᵢ_mean = map(duᵢ -> sum(integrate(t, -t .* duᵢ)), du)

        # proportion = tᵢ_mean[2] / (tᵢ_mean[1] + tᵢ_mean[2])



        # return tᵢ_mean, proportion
end

##
plot()
clr = [:blue :red]
mean_each_species, P_activated, plt = plot_mean_of_each_species(tmc_b_M_num,ks_tmc_b,tmc_b,TMC_B;species=["τ_weak" "τ_strong"])
plt
##
T2 = mean_each_species[:,2] ./ P_activated
plot!(0:40,T2,c=:green)
proportion_act = s_w_tmc_b[1:10:end]
mean_each_species[:,2] ./ proportion_act

function plot_P_act_v_t(M,k,data,binned;species=[],dt=0.01,dF=0.5)

        # ts = collect(0.1:0.1:1)
        ts = collect(0:0.01:3.0)


        # ts = 10 .^(range(-3,stop=1,length=5))
        # ts = [0,0.0001,0.001,0.01,0.1,1,2]
        Fs = 0:dF:maximum(data.F)
        # Fs = 0:dF:100

        M.solver = Tsit5  #change solver to Tsit5 (can handle ± std)
        sols = ftsolve(M,k,Fs,ts)
        ratios = []
        P2s = []
        for i in 1:length(sols)
                tmp_sol = sols[i]
                _u = tmp_sol(ts) 
                u = [_u[i,:] for i in 1:size(_u)[1]]
                ratio = u[2] ./ (u[1] .+ u[2])
                p2 = u[2]
                push!(ratios,ratio)
                push!(P2s,p2)
        end
        _ts = fill(collect(ts),length(Fs))
        
        # plot(ts,ratios,lw=3,color = :turbo, line_z = Fs',label="",colorbar_title ="F")
        # plot(ts,P2s,lw=3,color = :turbo, line_z = Fs',label="")
        # ratios_v_F = reduce(hcat,ratios)'

        # plot(Fs,reverse(ratios_v_F,dims=2),lw=3,color = :turbo, line_z = reverse(log10.(ts')),label="")
        # plot(Fs,ratios_v_F,lw=3,color = :turbo, line_z = log10.(ts'),label="",colorbar_title="log10(t)")

        P2s_v_F = reduce(hcat,P2s)
        @show maximum(P2s_v_F)
        # plot(Fs,P2s_v_F,lw=3,color = :turbo, line_z = log10.(ts'),label="",colorbar_title="log10(t)")

        # surface(Fs,ts,P2s_v_F',ylabel="t (s)",xlabel="Force (pN)",
        #                 colorbar_title ="P(strong)",#colorbar_titlefontsize=10,
        #                 fontsize=10)
        # SURFACE PLOT FOR Kaitao
        # Plots.scalefontsizes(2)
        lb = minimum(P2s_v_F)
        @show size(P2s_v_F)
        #append row of P = 0.5 at t > 3 to make sure that the colorbar goes from 0 to 0.5
        _P2s_v_F = vcat(P2s_v_F, ones(size(P2s_v_F)[2])'*0.5 )
        c1 = cgrad(:turbo, rev = false)
        p1 = plot(Fs,[ts;3.01],_P2s_v_F,st=:contour,
                        fill=(true),color=c1,
                        camera=(10,45),
                        clims=(0,0.5),
                        ylabel="Lifetime (s)",xlabel="Force (pN)",
                        colorbar_title ="P( (TCR-pMHC-CD4)* | t)",#colorbar_titlefontsize=10,
                        # colorbar_title ="P( (TCR-pMHC)* | t)",
                        # colorbar_title ="P( (j)* | t)",
                        fontsize=10)
        # xlabel!("F")
        # xlabel!("t")
        # ylabel!("P(strong)")
        # ylabel!("Proportion_strong = P(strong)/P(weak+strong)")
        ylims!((0,3))
        # Plots.scalefontsizes(1/2)
        # return p1
end
        ##
# pyplot()
# plotly()
# plot()


tmc_b_M_num = M(bimolecular_diss_bell!,u₀,u₀_type,tspan,[],ks_tm_b,0.0,Rodas5)
plot_mean_ft_fit(tmc_b_M_num,k_tmc_b,tmc_b,TMC_B;species="TCR-MHC-CD4")
tm_b_M_num = M(bimolecular_diss_bell!,u₀,u₀_type,tspan,[],[],0.0,Rodas5)
plot_mean_ft_fit(tm_b_M_num,k_tm_b,tm_b,TM_B;species="TCR-MHC")

Plots.scalefontsizes(2)
plot_P_act_v_t(tmc_b_M_num,ks_tmc_b,tmc_b,TMC_B;species=["τ_weak" "τ_strong"])
# savefig(figures_dir*"P_activated_TMC.png"); savefig(figures_dir*"P_activated_TMC.svg")

plot_P_act_v_t(tm_b_M_num,ks_tm_b,tm_b,TM_B;species=["τ_weak" "τ_strong"])
# savefig(figures_dir*"P_activated_TM.png"); savefig(figures_dir*"P_activated_TM.svg")


function plot_Proportion_act_v_t(M,k,data,binned;species=[],dt=0.01,dF=0.5)

        ts = collect(0.:0.01:1)
        # ts = collect(0:0.01:3.0)


        # ts = 10 .^(range(-3,stop=1,length=5))
        # ts = [0,0.0001,0.001,0.01,0.1,1,2]
        Fs = 0:3:maximum(41)
        # Fs = 0:dF:100

        M.solver = Tsit5  #change solver to Tsit5 (can handle ± std)
        sols = ftsolve(M,k,Fs,ts)
        ratios = []
        P2s = []
        for i in 1:length(sols)
                tmp_sol = sols[i]
                _u = tmp_sol(ts) 
                u = [_u[i,:] for i in 1:size(_u)[1]]
                ratio = u[2] ./ (u[1] .+ u[2])
                p2 = u[2]
                push!(ratios,ratio)
                push!(P2s,p2)
        end
        _ts = fill(collect(ts),length(Fs))
        
        plt= plot(ts,ratios,lw=3,color = :turbo, line_z = Fs',
                                # xaxis=:log,  
                                clims=(0,40),                      
                                label="",
                                fontsize=10,
                                colorbar_title ="Force (pN)",colorbar_titlefontsize=20)
        # xlims!((0,1))
        xlabel!("Lifetime (s)")
        # ylabel!("P((TCR-pMHC-CD4)*) / ΣP")
        ylabel!("P((TCR-pMHC)*) / ΣP")

        # plot(ts,P2s,lw=3,color = :turbo, line_z = Fs',label="")
        ratios_v_F = reduce(hcat,ratios)
        @show size(ratios_v_F)
        # plot(Fs,reverse(ratios_v_F,dims=2),lw=3,color = :turbo, line_z = reverse(log10.(ts')),label="")
        # plot(Fs,ratios_v_F,lw=3,color = :turbo, line_z = log10.(ts'),label="",colorbar_title="log10(t)")

        P2s_v_F = reduce(hcat,P2s)
        # @show maximum(P2s_v_F)
        # plot(Fs,ratios_v_F',lw=3,color = :turbo, line_z = log10.(ts'),label="",colorbar_title="log10(t)")
        # plot(0:10)
        # c1 = cgrad(:turbo, rev = false)

        # plot(Fs,ts,reverse(ratios_v_F,dims=1),st=:contour,
        #                 fill=(true),color=c1)#,
                        # ylabel="t (s)",xlabel="Force (pN)",
                        # colorbar_title ="P(strong) / (P(strong) + P(weak))",#colorbar_titlefontsize=10,
                        # fontsize=10)
        # SURFACE PLOT FOR Kaitao
        # Plots.scalefontsizes(2)
        # lb = minimum(P2s_v_F)
        # @show size(P2s_v_F)
        #append row of P = 0.5 at t > 3 to make sure that the colorbar goes from 0 to 0.5
        # _P2s_v_F = vcat(P2s_v_F, ones(size(P2s_v_F)[2])'*0.5 )
        # c1 = cgrad(:turbo, rev = false)
        # p1 = plot(Fs,ts,ratios_v_F,st=:contour,
        #                 fill=(true),color=c1,
        #                 camera=(10,45),
        #                 clims=(0,1),
        #                 ylabel="t (s)",xlabel="Force (pN)",
        #                 colorbar_title ="P( (j)* | t)",#colorbar_titlefontsize=10,
        #                 fontsize=10)
        # xlabel!("F")
        # xlabel!("t")
        # ylabel!("P(strong)")
        # ylabel!("Proportion_strong = P(strong)/P(weak+strong)")
        # ylims!((0,3))
        # Plots.scalefontsizes(1/2)
        return plt, ratios
end

pyplot()
Plots.scalefontsizes(1/2)
Plots.scalefontsizes(1.8)
plt_tm, P_activated_v_F_and_t_tm = plot_Proportion_act_v_t(tm_b_M_num,ks_tm_b,tm_b,TM_B;species=["τ_weak" "τ_strong"])
plt_tm
# savefig(figures_dir*"Proportion_activated_TM.png"); savefig(figures_dir*"Proportion_activated_TM.svg")

plt_tmc, P_activated_v_F_and_t_tmc = plot_Proportion_act_v_t(tmc_b_M_num,ks_tmc_b,tmc_b,TMC_B;species=["τ_weak" "τ_strong"])
plt_tmc
# savefig(figures_dir*"Proportion_activated_TMC.png"); savefig(figures_dir*"Proportion_activated_TMC.svg")



P_activated_v_F_and_t_tm = reduce(hcat,P_activated_v_F_and_t_tm)
P_activated_v_F_and_t_tmc = reduce(hcat,P_activated_v_F_and_t_tmc)

CSV.write(res_dir*"P_activated_v_F_and_t_tm.csv",P_activated_v_F_and_t_tm)
CSV.write(res_dir*"P_activated_v_F_and_t_tmc.csv",P_activated_v_F_and_t_tmc) 

using Tables

CSV.write(res_dir*"P_activated_v_F_and_t_tm.csv", 
                DataFrame(P_activated_v_F_and_t_tm, :auto))

CSV.write(res_dir*"P_activated_v_F_and_t_tmc.csv", 
                DataFrame(P_activated_v_F_and_t_tmc, :auto))
