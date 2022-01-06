using Optim: initial_state

# module model_fit

using DifferentialEquations, Flux, Optim, DiffEqFlux, DiffEqSensitivity, Plots
using Zygote, ForwardDiff, ReverseDiff
using BlackBoxOptim#, GalacticOptim
using BenchmarkTools
using RecursiveArrayTools
# -------------------------------------------------------------------------


mutable struct M_global
    model   # model function
    u₀      # initial conditions
    tspan   # time span
    dens    # site densities
    cons    # constant parameters
end

mutable struct DissData
    y
    t
    f
end


function run_model(M,p,dens,cons,t) # version of run function with multiple models
  # f = (du,u,p,t) -> model(du,u,10 .^ p,t,densities,constants)
  f = (du,u,p,t) -> M.model(du,u,p,t,dens,cons)
  tspan = (0.0,maximum(t))
  # # tmp_prob = ODEProblem(f,M.u₀,M.tspan,p) #allows to specify tspan
  # tmp_prob = ODEProblem(f,M.u₀,tspan,p)  #automatic tspan

  # # uncomment to fit the initial conditions as well...
  # # #u₀ = [(1-p[1])/2, (1-p[1])/2, p[1]]
  # u₀ = [(1-p[1])*0.610, (1-p[1])*0.390, p[1]] #ratio depends on γ and δ affinities
  # p = p[2:end]
  u₀, p = initial_u(p,M) 
  tmp_prob = ODEProblem(f,u₀,tspan,p)
  # # end fit initial conditions
  tmp_sol = solve(tmp_prob,AutoVern7(KenCarp4()),saveat=t, abstol=1e-8,reltol=1e-8,verbose=false)
end

function loss(p,data,M,unique_t)
  Δy_all = [] # keep track of errors for all curves (change this to array of arrays!)
  # Δy_all = Zygote.Buffer([[0.,0.],[0.,0.]],length(data))
  for i in 1:length(data.y) # for each bin
    # run model with densities and parameter/dens/cons combinations
    tmp_sol = run_model(M, p, M.dens[i], M.cons[i], unique_t[i])

    if tmp_sol.retcode != :Success #if there is an error return infitite loss
      Δy_sol = fill(Inf, length(data.t[i]))
    else
      Δy_sol = find_diffs(tmp_sol,data.t[i],unique_t[i],data.y[i])
    end
    Δy_all = vcat(Δy_all,Δy_sol)
  end

  t_all = reduce(vcat,data.t)

  Losses = Δy_all .^2 #.* t_all #.^2
  L = sum(Losses) # calculate error across all differences from all curves

  return L
end


function global_LS_fit(M,data,options)
# fit function - given model, data, densities, options get best parameters
    unique_t = [unique(times) for times in data.t]
    loss_in = (p) -> loss(p,data,M,unique_t)
    callback_in = (p,l) -> callback_global(p,l,data,M,unique_t)

    # lb,ub = options
    bounds = options

    @show loss_in(pinit)

    # sr = [(0.0,1.0),(0.0,100.0)]
    # res_bbo = bboptimize(loss_in; SearchRange = sr,
    # res_bbo = bboptimize(loss_in; SearchRange = options,
    res_bbo = bboptimize(loss_in; SearchRange = bounds,
                                NumDimensions = length(pinit),
                                Method = :adaptive_de_rand_1_bin_radiuslimited,
                                # Method = :separable_nes,
                                NThreads=Threads.nthreads()-1,
                                MaxSteps = 20000)#,

    p₀_bbo = best_candidate(res_bbo)
    @show p₀_bbo
    # p₀_bbo =[ 5.045339922925327, 0.2110701531504992]
    lb,ub = boundaries(bounds)
    res = optimize(loss_in,lb,ub,p₀_bbo, Fminbox(BFGS()), autodiff = :forward,
                                Optim.Options(show_trace=true,
                                              f_tol = 1e-8,
                                              outer_iterations = 10))

    rates = res.minimizer
    l = res.minimum
    # bootstrap_rates = bootstrapping(rates,data,model,dens,cons,t
    hes = ForwardDiff.hessian(loss_in,rates)

    println("Model used == ", M.model)
    println("Loss (MSE) == ", l)
    println("Rates == ", rates)
    res = OptRes(rates, l, hes)
    return res
end



function find_diffs(tmp_sol,t,unique_t,y) # solutions use unique time points,
# need to evaluate the solution at each experimental time point as well as
# calculating the difference between fit and experiment

    Σ_sol = sum(Array(tmp_sol),dims=1)

    y = log.(y) # remove for adhesion
    Σ_sol = log.(abs.(Σ_sol)) # remove for adhesion


    Δy = []
    for (i, tᵢ) in enumerate(t)
        ind = isequal.(unique_t,tᵢ)
        yᵢ_exp = Σ_sol[ind][1]
        yᵢ_obs = y[i]
        Δyᵢ = yᵢ_obs .- yᵢ_exp
        Δy=vcat(Δy,Δyᵢ)
    end
    return copy(Δy) #Δy[1:end] # copy(Δy)
end

# function callback(p,l,tmp_sol,data,model,dens,cons,t)
function callback_global(p,l,data,M,unique_t)
  # @show l, p
  @show l

  false
end

function initial_u(k,M::M_global)
  if M.model == cd3_dissociation!
    # need to pass first parameter, k[1], as initial trimolecular population.
    # in this problem , we don't know the initial trimolecular population, hence 
    # we can keep it as a parameter that can be fitted.
    u₀ = [(1-k[1])*u₀_bi[1], (1-k[1])*u₀_bi[2], k[1]] #ratio depends on γ and δ affinities
    
    # all other paramters, aside from k[1], are chemical rates, so return them as rates
    k = k[2:end]
  elseif M.model == trimolecular_diss!
    u₀ = [1-k[3], k[3]]
    k = k[1:2]
  elseif M.model == trimolecular_diss_no_on2!
    u₀ = [1-k[2], k[2]]
    k = k[1]
  elseif M.model == bimolecular_diss_two_states!
    u₀ = [1-k[3], k[3]]
    k = k[1:2]
  else
    u₀ = M.u₀
  end
  u₀ = convert.(eltype(k),u₀)

  return u₀, k
end

hessian2σ_LS(H,l,M,N) = sqrt.(abs.( diag(inv(H)) )) * l #/(M-N)
