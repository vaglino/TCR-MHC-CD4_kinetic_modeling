using ColorSchemes: length
using ColorSchemes

# solve model for given rates and force interval
function ftsolve(M,rates,Fs,ts)
    sols = []
    for f in Fs
        # M.cons = f
        M.f = f
        sol = solve_prob(M,rates;tᵢ=ts)
        push!(sols,sol)
    end
    return sols
end

# take model solution and construct P and PDF surfaces for plotting
function construct_surf(sols,Fs,ts)
    p = []
    pdf = []
    for (i,f) in enumerate(Fs)

        pᵢ = sum(Array(sols[i](ts))',dims=2)
        pdfᵢ = - sum( sols[i](ts, Val{1}) , dims=2) 
        
        p = push!(p,pᵢ)
        pdf = push!(pdf,pdfᵢ)
    end
    p = reduce(hcat,p)
    pdf = reduce(hcat,pdf)
    return p, pdf
end

# compute pdf from p distribution
function p2pdf(sol)
    # mean_lifetime(sol)
    pdf = p ./ mean_lifetime(sol)
    pdf_ = sol(:,)
end


function plot_ft_fit(M,rates,data;dt=0.1,dF=1.0)

    ts = 0:dt:maximum(data.t)
    Fs = 0:dF:maximum(data.F)
    # Fs = 0:dF:20
    sols = ftsolve(M,rates,Fs,ts)

    p, pdf = construct_surf(sols,Fs,ts)

    # plot contours

    p1 = plot(Fs,ts,log.(p),st=:contour,
                            fill=(true),#cgrad(:plasma)),
                            camera=(0,90),
                            clims=(-7,0),
                            ylabel="<t>",xlabel="F",clabel="ln(p)")
    scatter!(data.F,data.t,c=:white,markeralpha=0.7)
    display(p1)

    pex = surface(Fs,ts,(p),#,st=:contour,
                            fill=(true),#cgrad(:plasma)),
                            # camera=(0,90),
                            # clims=(-7,0),
                            zlabel="p",xlabel="F")
    display(pex)

    p2 = plot(Fs,ts,p,st=:contour,
                        fill=(true),#,cgrad(:plasma)),
                        camera=(0,90),
                        ylabel="p",xlabel="F")
    display(p2)

    p3 = plot(Fs,ts,log.(pdf),st=:contour,
                        fill=(true,cgrad(:plasma)),
                        camera=(0,90),
                        ylabel="pdf",xlabel="F")

    # plot lnp curves

    p4=plot()
    lnp = p2lnp.(p)
    for i in 1:1:size(lnp,2)
            plot!(ts,lnp[:,i],c = i,
                                xlims=(0,ts[end]),
                                ylims=(-7,0),
                                legend=false,
                                ylabel="ln(p)",xlabel="t")
    end


    # plot mean lifetime distribution
    t_means = mean_lifetime_from_S.(sols)
    p5 = plot(Fs,t_means,ylabel="<t>",xlabel="F")

    display(plot(p1,p2,p3,p4,p5))
    # display(p2)
    # display(p3)
    # display(p4)
    # display(p5)
    return sols
end


function plot_mean_ft_fit(M,k,data,binned;species=[],dt=0.1,dF=1.0)

    ts = 0:dt:M.tspan[2]
    Fs = 0:dF:maximum(data.F)

    M.solver = Tsit5  #change solver to Tsit5 (can handle ± std)
    sols = ftsolve(M,k,Fs,ts)
    mean_t = mean_lifetime_from_S.(sols)
    plot!(Fs,mean_t,c=clr,
                    # ribbon=0.1,
                    markerstrokecolor=:auto,
                    linealpha = 1,
                    ylabel="<t>",
                    xlabel="F",
                    label=species)
    # plot_F_vs_t(clean_dissociation_data(binned))
    plot_F_vs_t(binned)
end

function extrapolate_fit(M,k,data,F_max;species=[],dt=0.1,dF=1.0)

    ts = 0:dt:M.tspan[2]
    Fs = 0:dF:maximum(F_max)

    # M.solver = Tsit5  #change solver to Tsit5 (can handle ± std)
    sols = ftsolve(M,k,Fs,ts)
    mean_t = mean_lifetime.(sols)
    h = plot!(Fs,log10.(mean_lifetime.(sols)),c=Color,
                    # ribbon=0.1,
                    markerstrokecolor=:auto,
                    linealpha = 1,
                    ylabel="<t>",
                    xlabel="F",
                    label=species)
end


function mean_ft_fit(M,k,data,F_max;species=[],dt=0.1,dF=1.0)

    ts = 0:dt:maximum(data.t)
    # ts = 0:dt:30
    Fs = 0:dF:F_max

    M.solver = Tsit5  #change solver to Tsit5 (can handle ± std)
    sols = ftsolve(M,k,Fs,ts)
    mean_t = mean_lifetime_linear.(sols)
end

##_____________________________________________________________________
# plotting for analytical model

function ftsolve_analytical(M,rates,Fs,ts)
    sols = []
    for f in Fs
        sol = map(t -> M.model(rates,M,t;f=f), ts)
        push!(sols,sol)
    end
    return sols
end
function construct_surf_analytical(sols,Fs,ts)
    p = []
    pdf = []
    for (i,f) in enumerate(Fs)

        pᵢ = sols[i] .* mean_lifetime_analytical(ts,sols[i])
        # @show mean_lifetime_analytical(sols[i])
        pdfᵢ = sols[i] 

        p = push!(p,pᵢ)
        pdf = push!(pdf,pdfᵢ)
    end
    p = reduce(hcat,p)
    pdf = reduce(hcat,pdf)
    return p, pdf
end

function mean_lifetime_analytical(t,M,k)
    # NOTE: gives inaccurate results when t_max is too short. Integrating the 
    # survival distribution instead of the pdf tends to lead to more
    # accurate results.
    # mean = integrate(t,t.*pdf)
    S = pdf2S_quad(M,k,t)
    mean = mean_lifetime_from_S(t,S)
    # @show mean
    return mean
end
function mean_lifetime_analytical(t,pdf::Array)
    # NOTE: gives inaccurate results when t_max is too short. Integrating the 
    # survival distribution instead of the pdf tends to lead to more
    # accurate results.
    mean = integrate(t,t.*pdf)
    # S = pdf2S(t,pdf)
    # mean = mean_lifetime_from_S(t,S)
    # @show mean
    return mean
end
function mean_ft_fit_analytical(M,k,data;F_max=30.,t_max=100.,dt=0.01,dF=1.0)

    t = 0:dt:t_max
    Fs = 0:dF:F_max

    # sols = ftsolve_analytical(M,k,Fs,ts)
    # # mean_t = map(sol -> mean_lifetime_analytical(ts,sol), sols)
    mean_t = []
    # for F in Fs
    #     M.f = F
    #     push!(mean_t, mean_lifetime_analytical(ts,M))
    # end
    for F in Fs
        # push!(mean_t,quadgk(t -> t*M.model(k,M,t;f=F), 0., t_max)[1])
        M.f = F
        push!(mean_t,mean_lifetime_analytical(t,M,k))
    end
    return mean_t
end

function mean_ft_fit_analytical_2(M,k,data;F_max=30.,t_max=100.,dt=0.01,dF=1.0)

    t = 0:dt:t_max
    Fs = 0:dF:F_max

    mean_t = []
    for F in Fs
        M.f = F
        push!(mean_t,quadgk(t -> t*M.model(k,M,t;f=F), 0., t_max)[1])
    end
    return mean_t
end

function survival_dist_analytical(M,k;F=0.,t_max=10.,dt=0.1)
    M.f = F
    t = 0:dt:t_max
    p = pdf2S_quad(M,k,t)
end


function survival_dist_2d_analytical(M,k;F_max=30.,t_max=10.,dF=0.5,dt=0.1)
    p =[]
    for F in 0:dF:F_max
        pᵢ = survival_dist_analytical(M,k;F=F,t_max=t_max,dt=dt)
        push!(p,pᵢ)
    end
    p = reduce(hcat,p)
end

function plot_survival_surface_analytical(M,k,data;F_max=30,t_max=10,dF=0.5,dt=0.1)
    p = survival_dist_2d_analytical(M,k;F_max=F_max,t_max=t_max,dF=dF,dt=dt)
    ts = 0:dt:t_max
    Fs = 0:dF:F_max

    lb = -10
    p[p .< exp(lb)] .= exp(lb)
    # Plots.scalefontsizes(2)

    c1 = cgrad(:turbo, rev = false)
    p1 = plot(Fs,ts,log.(p),st=:contour,
                            fill=(true),color=c1,#cgrad(:plasma)),
                            camera=(0,90),
                            clims=(lb,0),
                            ylabel="Lifetime  (s)",xlabel="Force (pN)",
                            colorbar_title ="ln(p)",colorbar_titlefontsize=12,
                            fontsize=10)
    scatter!(data[1],data[2],c=:white,
                markersize=3,markeralpha=0.7,markerstrokealpha=0.3,
                label=nothing,thickness_scaling = 2, grid=false)
    # Plots.scalefontsizes(1/2)
    return p1
end


# figures_dir = "C:/Users/stravaglino3/Downloads/julia_code/CD3gd/figures/"

# pyplot()
# plot()
# p1 = plot_survival_surface_analytical(γϵ_ft_M,ks_γϵ,γϵ_no_boot;F_max=30,t_max=15,dF=0.5,dt=0.1)
# savefig(figures_dir*"survival_surface_γϵ.png")

# plot()
# p2 = plot_survival_surface_analytical(δϵ_ft_M,ks_δϵ,γϵ_no_boot;F_max=30,t_max=15,dF=0.5,dt=0.1)
# savefig(figures_dir*"survival_surface_δϵ.png")

# plot()
# p3 = plot_survival_surface_analytical(mix_ft_M,ks_mix,mix_no_boot;F_max=30,t_max=60,dF=0.5,dt=0.1)
# savefig(figures_dir*"survival_surface_mix.png")




function trimolecular_proportion_vs_t(M_bi1,M_bi2,M_tri,k1,k2,k3;F=0,t_max=10)
    p_bi1 = survival_dist_analytical(M_bi1,k1; F=F, t_max=t_max,dt=0.05) * u₀_mix[1]
    p_bi2 = survival_dist_analytical(M_bi2,k2; F=F, t_max=t_max,dt=0.05) * u₀_mix[2]
    p_tri = survival_dist_analytical(M_tri,k3; F=F, t_max=t_max,dt=0.05) * u₀_mix[3]

    ratio = p_tri ./ ((p_bi1+p_bi2) + p_tri)
end

# trimolecular_proportion_vs_t(γϵ_ft_M,δϵ_ft_M,mix_ft_M,ks_γϵ,ks_δϵ,ks_mix)

function trimolecular_proportion_v_F_and_t(M_bi1,M_bi2,M_tri,k1,k2,k3;Fs=0:30,t_max=10)
    ratios = []
    for F in Fs
        ratio_at_F = trimolecular_proportion_vs_t(M_bi1,M_bi2,M_tri,k1,k2,k3;F=F,t_max=10)
        push!(ratios,ratio_at_F)
    end
    return ratios
end


# ratios = trimolecular_proportion_v_F_and_t(γϵ_ft_M_w_uncert,δϵ_ft_M,mix_ft_M,ks_γϵ,ks_δϵ,ks_mix;Fs=Fs)
# CSV.write(results_dir*"tri_proportion/0pN.csv", measurements2Df(0:0.05:10,ratios[1];labels=["p_tri","std"]))
# CSV.write(results_dir*"tri_proportion/7pN.csv", measurements2Df(0:0.05:10,ratios[2];labels=["p_tri","std"]))
# CSV.write(results_dir*"tri_proportion/13pN.csv", measurements2Df(0:0.05:10,ratios[3];labels=["p_tri","std"]))
# CSV.write(results_dir*"tri_proportion/20pN.csv", measurements2Df(0:0.05:10,ratios[4];labels=["p_tri","std"]))

# plot(0:0.05:10,ratios,lw=3,color = :turbo, line_z = (1:length(ratios))')
# plot(0:0.05:10,ratios[1],lw=3)
# plot!(0:0.05:10,ratios[2],lw=3)
# plot!(0:0.05:10,ratios[3],lw=3)
# plot!(0:0.05:10,ratios[4],lw=3)

# ylabel!("P(Tri)")
# xlabel!("Bond Lifetime (s)")

function plot_survival_fit(M,k,data; F=0, t_max=10, dt=0.1, clr=:blue)

    p = survival_dist_analytical(M,k;F=F,t_max=t_max,dt=dt)
    lnp = log.(abs.(p))

    # =
    p_ = plot!(0:dt:t_max,lnp,color=clr,lw=5)
    scatter!(data[1],data[2],color=clr)
    # =#

    # lnn = normalize_lnp_fit(lnp,data[1])
    # p_ = plot!(0:dt:t_max,lnn,color=clr)
    # lnn_data = normalize_lnp(data[2])
    # scatter!(data[1],lnn_data,color=clr)
    # lnp = lnn

    return lnp, p_
    # end
end


function survival_fit_v_F(M,k,data; F=0, t_max=10, dt=0.1, clr=:blue)
    bins, binned_t, edges = bin_dissociations(data, 2; bin_type = "force", edges=[0,F-4,F+4])
    bin_data = [sort(binned_t[2]),calculate_lnp(binned_t[2])]

    lnp, p_ = plot_survival_fit(M,k,bin_data; F=F, t_max=t_max, dt=dt, clr=clr)
    plot!(p_)
    return lnp,bin_data,p_
end


function survival_multiple_bins(M,k,data; Fs=0:5:30, t_max=10, dt=0.1, clr=:blue)
    lnps=[]
    bins=[]
    p1 = plot()
    clr = palette(:turbo,4)
    for (i,F) in enumerate(Fs)
        lnp,bin_data,p_ = survival_fit_v_F(M,k,data; F=F, t_max=t_max, dt=dt, clr=clr[i])
        plot!(p_)
        push!(lnps,lnp)
        push!(bins,bin_data)
    end
    return lnps,bins,p1
end

function survival2Df(t,lnps,Fs)
    df = DataFrame()
    df[!, "t (s)"] = t
    for (lnp,F) in zip(lnps,Fs)
        fit = Measurements.value.(lnp)
        std = Measurements.uncertainty.(lnp)
        label = "lnp "*string(F)* "pN"
        label_std = "std "*string(F)* "pN"
        # df.F = lnp
        # df."σ" = std 
        df[!, label] = fit
        df[!, label_std] = std
        # push!(df,lnp)
        # push!(df,σ )
    end
    return df
end

# gr()
# # plotly()
# Fs = [0,7.5,15,22.5]
# Fs = [0,7,13,20]

# plot()
# lnps_γϵ,bins_γϵ,p_γϵ_ = survival_multiple_bins(γϵ_ft_M,k_γϵ,γϵ_no_boot; Fs=Fs, t_max=10, dt=0.1, clr=:blue) # [ΓΕ.F, ΓΕ.t]
# plot(p_γϵ_)
# ylims!((-6,0))
# ylims!((0,1))
# plot()
# lnps_δϵ,bins_δϵ,p_δϵ_ = survival_multiple_bins(δϵ_ft_M,k_δϵ,δϵ_no_boot; Fs=Fs, t_max=10, dt=0.1, clr=:blue)  # [ΔΕ.F, ΔΕ.t]
# plot(p_δϵ_)
# ylims!((-6,0))
# plot()
# lnps_mix,bins_mix,p_mix_ = survival_multiple_bins(mix_ft_M,k_mix,mix_no_boot; Fs=Fs, t_max=60, dt=0.1, clr=:blue) # [Mix.F, Mix.t]
# plot(p_mix_)
# ylims!((-6,0))


# df_lnp_γϵ = survival2Df(0:0.1:10,lnps_γϵ,Fs)
# CSV.write(res_dir*"lnp_γϵ.csv", df_lnp_γϵ)
# df_lnp_δϵ = survival2Df(0:0.1:10,lnps_δϵ,Fs)
# CSV.write(res_dir*"lnp_δϵ.csv", df_lnp_δϵ)
# df_lnp_mix= survival2Df(0:0.1:60,lnps_mix,Fs)
# CSV.write(res_dir*"lnp_mix.csv", df_lnp_mix)



function survival_dist_2Df(bins,Fs)
    df = DataFrame()
    n_max = maximum( map(x -> length(x[1]), bins) )

    for (bin,F) in zip(bins,Fs)
        
        n_bin = length(bin[1])

        label_t = "t (s) at "*string(F)* "pN"
        label_lnp = "ln(p) at "*string(F)* "pN"
        # add columns to df and pad so they are all the same length
        df[!, label_t] =  [bin[1]; [missing for _ in 1:(n_max-n_bin)]]
        df[!, label_lnp] = [bin[2]; [missing for _ in 1:(n_max-n_bin)]]
    end
    return df
end

# df_survival_γϵ = survival_dist_2Df(bins_γϵ,Fs)
# df_survival_δϵ = survival_dist_2Df(bins_δϵ,Fs)
# df_survival_mix = survival_dist_2Df(bins_mix,Fs)
# CSV.write(results_dir*"survival_dist_γϵ.csv", df_survival_γϵ)
# CSV.write(results_dir*"survival_dist_δϵ.csv", df_survival_δϵ)
# CSV.write(results_dir*"survival_dist_mix.csv", df_survival_mix)


