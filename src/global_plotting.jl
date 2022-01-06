function plot_sum(rates,Color=nothing)
    prob = ODEProblem(f,u₀,tspan,rates)
    # prob = remake(prob, p=10 .^ rates)
    sol = solve(prob,Tsit5(),saveat=0.1)
    summ = reduce(vcat,sum(sol,dims=1))
    if Color == nothing
        plot!(sol.t,summ)
    else
        display(plot!(sol.t,summ,c=Color))
    end
    return DataFrame(t = sol.t,n = summ)
end

function plot_n_fit(_M,p,dens,cons;Color=nothing)

    M = deepcopy(_M)
    M.solver = Tsit5
    t = M.tspan[1]:0.1:M.tspan[2]

    if first(methods(M.model)).nargs == 5
        tmp_sol = map(tᵢ -> M.model(p,tᵢ,dens,cons), t)'
        Σ_sol = reduce(vcat,sum(tmp_sol,dims=1))
    else
        tmp_sol = run_model(M,p,dens,cons,t)
        Σ_sol = reduce(vcat,sum(tmp_sol,dims=1)) ./ (dens[1]*dens[2])
    end

    Pa = n2Pa( Σ_sol * (dens[1]*dens[2]) ) 

    if Color == nothing
        plot!(t,Σ_sol)
    else
        display(plot!(t,Σ_sol,c=Color,lw=5))
    end
    # @show typeof(Σ_sol)
    if typeof(Σ_sol[1]) <: Measurement
        df_n = measurements2Df(t,Σ_sol;labels=["n","n_std"])
        df_Pa = measurements2Df(t,Pa;labels=["Pa","Pa_std"])
        df = innerjoin(df_n, df_Pa, on = :t)

    else
        df_n = DataFrame(t = t, n = Σ_sol)
        df_Pa = DataFrame(t = t, Pa = Pa)
        df = innerjoin(df_n, df_Pa, on = :t)
    end

    return df
end


function display_rates(rates)
    for rate in rates
        @show rate
    end
end
