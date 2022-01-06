using DifferentialEquations, Optim, DiffEqFlux, DiffEqSensitivity, Plots
using Zygote, ForwardDiff, ReverseDiff
# using BlackBoxOptim#,GalacticOptim
using BenchmarkTools
using RecursiveArrayTools
# -------------------------------------------------------------------------

mutable struct F_t_Data
    F # forces vector
    t # lifetimes vector
    kd # kernel density of F
end

mutable struct M
    model   # model function
    u₀      # initial conditions
    u₀_type # string defining initial condtion type
    tspan   # time span
    dens    # site densities
    cons    # constant parameters
    f       # force
    solver  # ode solver
end

mutable struct OptRes
    rates
    l
    hes
end

function solve_prob(M,p;tᵢ=[])

    f = (du,u,p,t) -> M.model(du,u,p,t,M.dens,M.cons;f=M.f)
    u₀ = initial_u(p,M)
    # @show u₀
    # u₀ = [p[end], 1 - p[end]]
    # p = p[1:end-1]

    prob = ODEProblem(f,u₀,M.tspan,p)
    # sol = solve(prob,M.solver(),saveat=tᵢ,#save_start=false,
    #                         isoutofdomain=(u,p,t)->any(x->x<0,u),
    #                         abstol=1e-10,
    #                         reltol=1e-10,
    #                         verbose=false)
    sol1 = solve(prob,M.solver(),tstops=tᵢ,dtmax=0.1,
                            # callback=PositiveDomain(),
                            isoutofdomain=(u,p,t)->any(x->x<0,u),
                            abstol=1e-10,
                            reltol=1e-10,
                            verbose=false)
    # h = plot(sol)
    # plot!(sol1)
    # display(h)
    return sol1
end

function mle_loss(p,data,M)

    L = []
    for (Fᵢ,tᵢ) in zip(data.F,data.t)
    # Threads.@threads for i in 1:length(data.F)

    #     Fᵢ = data.F[i]
    #     tᵢ = data.t[i]

        M.f = Fᵢ
        # @show Fᵢ, tᵢ
        # check whether model has analytical solution
        # if M.model == two_state_catch_bond_analytical_w_cons || M.model == two_state_catch_bond_analytical 
        # @show first(methods(M.model)).nargs
        if first(methods(M.model)).nargs == 4
            Lᵢ = M.model(p,M,tᵢ;f=Fᵢ)
            push!(L,Lᵢ)
            if Lᵢ < 0.0
                @show Fᵢ,tᵢ
            end

        else # solve ode system
            sol = solve_prob(M,p;tᵢ)
            #check whether integration is successful
            if sol.retcode != :Success || sol.t[end] < tᵢ
                push!(L,0.0) # Likelihood of model is 0
            else
                # ŷᵢ = sum(sol(tᵢ))

                # the Likelihood (pdf) of a survival distribution S is -dS/dt 
                # extract derivative of S, by querying it at time tᵢ with Val{1}
                # then use linear property of derirvative dΣᵢ(nᵢ) = dn₁+...+dnₙ
                Lᵢ = - sum( sol(tᵢ, Val{1}) ) 
                if Lᵢ < 0
                    # @show Lᵢ, mean_t, ŷᵢ, tᵢ
                end
                push!(L,Lᵢ)
            end
        end
    end
    logL = logLikelihood(L)
    # try
        # logL = wLogLikelihood(L,1 ./ pdf(data.kd,data.F))
    # catch
    #     @show L
    # end
    # # sum(L .< 0.) == 0 || println(sum(L .< 0.) )
    # logL = wLogLikelihood(L,1 ./ pdf(data.kd,data.F))
    # logL = wLogLikelihood(L,data.t .^ (0.05) )
    # logL = wLogLikelihood(L,log.(data.t .+ 1.0) ./ pdf(data.kd,data.F))
    # logL = wLogLikelihood(L,data.t .^(0.1) ./ pdf(data.kd,data.F))

    # CSV.write("L_an.csv",DataFrame(L = L))
    return logL
end

# log likelihood
function logLikelihood(L)
    logL = -sum(log.(L))
end

# weighted log likelihood
function wLogLikelihood(L,w)
    wlogL = -sum(w .* log.(L))
end

# calculates mean lifetime from the survival distrubution S (or p) which is provided as DiffEq solution struct
using NumericalIntegration
function mean_lifetime_from_S(sol)
    t = sol.t
    u = sum(sol,dims=1)'
    sum(integrate(t,u))
end
function mean_lifetime_from_S(t,u)
    integrate(t,u)
end


pdf2S(t,pdf) = 1 .- cumul_integrate(t, pdf)
S2pdf(t,p) = -diff(p)./diff(t)
S2pdf(S) = - vec(sum( sol(S.t, Val{1}), dims=1))

function pdf2S_quad(M,k,t)
    F_and_errors = quadgk.(t -> M.model(k,M,t;f=M.f), 0., t)
    F = map(x->x[1],F_and_errors)
    S = 1 .- F
end

#________________________________________________________________________________________________
# Fitting routine
# first, explore parameter space with differential evolution algorithm, in log-space,
# then use best fit rates as initial position for further BFGS optimization to hone in the 
# minimum, in linear space.
function mle_fit(M,data,options)
    bounds = options

    loss_in = (p) -> mle_loss(10 .^ p,data,M)

    loss_in(pinit)

    # BBO optimization (use log bounds for easier search across orders of magnitude)
    # #= -----------------------------------------------------------------------
    log_bounds = map(x-> log10.( (x[1].+1e-7, x[2]) ) , bounds)
    res_bbo = bboptimize(loss_in; SearchRange = log_bounds,
    # res_bbo = bboptimize(loss_in; SearchRange = bounds,
                                NumDimensions = length(pinit),
                                Method = :adaptive_de_rand_1_bin,
                                # Method = :adaptive_de_rand_1_bin_radiuslimited,
                                # Method = :separable_nes,
                                NThreads=Threads.nthreads()-1,
                                MaxSteps = 5000)#,
                                # MaxSteps = 1500)#,
    #
    p₀_bbo = best_candidate(res_bbo)
    p₀_bbo = 10 .^ p₀_bbo

    # p₀_bbo = pinit

    # then hone in on minimum with gradient based method
    loss_in = (p) -> mle_loss(p,data,M)

    # bounds = map(x->log10.(x), options)
    lb,ub = boundaries(bounds)

    od = OnceDifferentiable(loss_in, p₀_bbo; autodiff = :forward);
    res = optimize(od,lb,ub,p₀_bbo,Fminbox(BFGS()),#,
    # res = optimize(od,p₀_bbo,BFGS(),
    # res = optimize(loss_in,lb,ub,p₀_bbo,Fminbox(NelderMead()),#
                                        Optim.Options(show_trace=true,
                                        iterations = 100,
                                        # iterations = 50,
                                        outer_iterations = 2))
                                        # f_tol = 1e-5,

    rates = res.minimizer
    l = res.minimum
    hes = ForwardDiff.hessian(loss_in,rates)

    println("Model used == ", M.model)
    println("Loss (MLE) == ", l)
    println("Rates == ", rates)
    res = OptRes(rates, l, hes)
    return res
end

callback = function (p, l)
  display(l)
  return false
end


using LinearAlgebra
# hessian2σ(H) = sqrt.(abs.( diag(inv(H)) ))
hessian2σ(H) = sqrt.( diag(inv(H)) )

hessian2σ(H,n) = sqrt.( diag(inv(H)) /n )

using QuadGK
function mean_lifetime(sol) # numerical integration to be changed
    # <t> = ∫p(t) between [0,∞), in this case <t> = ∫sol
    # stop inegration when solution reaches 0
    endpt_id = findfirst(iszero,sol[:])
    if endpt_id == nothing
        endpt = sol.t[end]
    else
        endpt = sol.t[endpt_id]
    end

    t_meanᵤ, err = quadgk(sol,0.0,endpt,atol=1e-10,rtol=1e-10)
    # if sum(t_meanᵤ) == 0.0
    #     @show sol
    # end
    t_mean = sum(t_meanᵤ)
end


#________________________________________________________________________________________________
# OLD optimization code, contains other ways to optimize, such as using GalacticOptim
#= 
function optimization(M,data,options)
    bounds = options

    loss_in = (p) -> mle_loss(10 .^ p,data,M)

    loss_in(pinit)

    # BBO optimization (use log bounds for easier search across orders of magnitude)
    # #= -----------------------------------------------------------------------
    log_bounds = map(x->log10.(x .+ 1e-5), bounds)
    res_bbo = bboptimize(loss_in; SearchRange = log_bounds,
    # res_bbo = bboptimize(loss_in; SearchRange = bounds,
                                NumDimensions = length(pinit),
                                Method = :adaptive_de_rand_1_bin,
                                # Method = :adaptive_de_rand_1_bin_radiuslimited,
                                # Method = :separable_nes,
                                NThreads=Threads.nthreads()-1,
                                MaxSteps = 1500)#,
    #
    p₀_bbo = best_candidate(res_bbo)
    # p₀_bbo = [0.809846, -4.85167, -0.371314, -1.18161, -1.02823, -3.21575, 1.12937, 0.366384]
    p₀_bbo = 10 .^ p₀_bbo
    #-----------------------------------------------------------------------
    # =#

    # @show p₀_bbo
    # p₀_bbo = log10.(pinit)
    # p₀_bbo = pinit

    #------------------------------------------------------------------------
    #------------------------------------------------------------------------


    # then hone in on minimum with gradient based method
    loss_in = (p) -> mle_loss(p,data,M)

    # bounds = map(x->log10.(x), options)
    lb,ub = boundaries(bounds)

    # @show so.minimize(x ->loss_in(x), pinit,method="BFGS")
    # @show so.minimize(loss_in, pinit,method="BFGS")
    # @show so.newton(loss_in, pinit)#,method="Newton")
    # #=
    od = OnceDifferentiable(loss_in, p₀_bbo; autodiff = :forward);
    res = optimize(od,lb,ub,p₀_bbo,Fminbox(BFGS()),#,
    # res = optimize(loss_in,lb,ub,p₀_bbo,Fminbox(NelderMead()),#
                                        Optim.Options(show_trace=true,
                                        iterations = 100,
                                        # iterations = 50,
                                        outer_iterations = 1))

                                        # f_tol = 1e-5,

    # =#
    # res = optimize(loss_in,p₀_bbo,lb,ub,Fminbox(BFGS()),
    #                                     iterations = 1,
    #                                     Optim.Options(iterations = 10))

    #= -----------------------------------------------------------------------
    optprob = GalacticOptim.OptimizationFunction((x,p) -> loss_in(x),p₀_bbo,GalacticOptim.AutoForwardDiff())
    # optprob = GalacticOptim.OptimizationFunction((x,p) -> loss_in(x))
    # # prob = GalacticOptim.OptimizationProblem(optprob,p₀_bbo)
    # # res = GalacticOptim.solve(prob,NelderMead(),cb = callback)
    #
    prob = GalacticOptim.OptimizationProblem(optprob,p₀_bbo;lb=lb,ub=ub)
    # )#
    res = GalacticOptim.solve(prob,BFGS(initial_stepnorm=0.01),
                                Fminbox(BFGS()),
                                cb = callback,
                                maxiters = 10)
    # res = GalacticOptim.solve(prob,BFGS(),cb = callback,maxiters = 20)
    =#

    rates = res.minimizer
    l = res.minimum

    hes = ForwardDiff.hessian(loss_in,rates)

    @show M.model
    @show l
    @show rates
    return rates, l, hes
end

=#