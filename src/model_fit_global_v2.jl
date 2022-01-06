using Revise, BenchmarkTools
using DifferentialEquations, Optim, DiffEqFlux, DiffEqSensitivity, Plots
using Zygote, ForwardDiff, ReverseDiff
using LinearAlgebra, RecursiveArrayTools, Measurements, Statistics
using DataFrames
# using BlackBoxOptim
# -------------------------------------------------------------------------

mutable struct Data
    n # number of bonds, or p
    t # lifetimes vector
end

mutable struct M_global
    model   # model function
    u₀      # initial conditions
    # u₀_type # string defining initial condtion type
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

function run_model(M,p,dens,cons,t) # version of run function with multiple models
    # f = (du,u,p,t) -> model(du,u,10 .^ p,t,densities,constants)
    f = (du,u,p,t) -> M.model(du,u,p,t,dens,cons)

    # u₀ = [1- p[3], p[3]]
    # p = p[1:2]
    # tmp_prob = ODEProblem(f,u₀,M.tspan,p)

    tmp_prob = ODEProblem(f,M.u₀,M.tspan,p)
    tmp_sol = solve(tmp_prob,M.solver(),saveat=t,
                                isoutofdomain=(u,p,t)->any(x->x<0,u),
                                abstol=1e-8,
                                reltol=1e-8,
                                verbose=false)
    # AutoVern7(KenCarp4()) works well
    # tmp_sol = solve(tmp_prob,Vern7(),saveat=t, abstol=1e-8,reltol=1e-8,verbose=false)
    tmp_sol
end


# function loss(p,data,M,unique_t)
#     Δy_all = [] # keep track of errors for all curves (change this to array of arrays!)
#     # Δy_all = Zygote.Buffer([[0.,0.],[0.,0.]],length(data))
#     for i in 1:length(data.n)
#       # run model with densities and parameter/dens/cons combinations
#       tmp_sol = run_model(M,p,M.dens[i],M.cons[i],unique_t[i])

#       # tmp_sol = run_delay_model(M,p,M.dens[i],M.cons[i],unique_t[i])

#       if tmp_sol.retcode != :Success #if there is an error return infitite loss
#         Δy_sol = fill(Inf, length(data.n[i]))
#       else
#         Δy_sol = find_diffs(tmp_sol,data.t[i],unique_t[i],data.n[i])
#       end
#       Δy_all = vcat(Δy_all,Δy_sol)
#       # Δy_all[i] = Δy_sol
#     end
#     # Δy = copy(Δy_all
#     # Losses = [sum(Δ .^ 2) for Δ in Δy]
#     # Losses = [sum( (Δy[i] .^ 2) .* t[i] .^2 ) for i in 1:length(Δy)]
#     t_all = reduce(vcat,data.t)

#     Losses = Δy_all .^2 .* t_all
#     # Losses = Δy_all .^2 .* t_all
#     L = sum(Losses) # calculate error across all differences from all curves

#     return L
# end

function loss(p,data,M,unique_t)
  Δy_all = [] # keep track of errors for all curves (change this to array of arrays!)
  # Δy_all = Zygote.Buffer([[0.,0.],[0.,0.]],length(data))
  for i in 1:length(data.n)
    # run model with densities and parameter/dens/cons combinations
    # if first(methods(M.model)).nargs == 4
    # if M.model == trimolecular_ana_zhu
    if first(methods(M.model)).nargs == 5
      tmp_sol = map(tᵢ -> M.model(p,tᵢ,M.dens[i],M.cons[i]), unique_t[i])'
      success = true
    else
      tmp_sol = run_model(M,p,M.dens[i],M.cons[i],unique_t[i])
      success = tmp_sol.retcode == :Success
      tmp_sol = Array(tmp_sol) ./ (M.dens[i][1] * M.dens[i][2])
    end

    # tmp_sol = run_delay_model(M,p,M.dens[i],M.cons[i],unique_t[i])

    if !success #if there is an error return infitite loss
      Δy_sol = fill(Inf, length(data.n[i]))
    else
      Δy_sol = find_diffs(tmp_sol,data.t[i],unique_t[i],data.n[i])
    end
    Δy_all = vcat(Δy_all,Δy_sol)
    # Δy_all[i] = Δy_sol
  end

  t_all = reduce(vcat,data.t)

  Losses = Δy_all .^2 #.* t_all
  # Losses = Δy_all .^2 .* t_all
  L = sum(Losses) # calculate error across all differences from all curves

  return L
end

function find_diffs(tmp_sol,t,unique_t,y) # solutions use unique time points,
# need to evaluate the solution at each experimental time point as well as
# calculating the difference between fit and experiment

    Σ_sol = sum(Array(tmp_sol),dims=1)

    # y = log.(y) # remove for adhesion
    # Σ_sol = log.(abs.(Σ_sol)) # remove for adhesion


    # Σ_sol = sum(tmp_sol,dims=1)
    # Δy = Zygote.Buffer(t,length(t))
    Δy = []
    # Δy = zeros(length(t))
    for (i, tᵢ) in enumerate(t)
        ind = isequal.(unique_t,tᵢ)
        yᵢ_exp = Σ_sol[ind][1]
        yᵢ_obs = y[i]
        Δyᵢ = yᵢ_obs .- yᵢ_exp
        # Δy[i] = Δyᵢ
        # push!(Δy,Δyᵢ)
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

function global_LS_fit(M,data,options)
  # fit function - given model, data, densities, options get best parameters
      unique_t = [unique(times) for times in data.t]
      loss_in = (p) -> loss(p,data,M,unique_t)
      callback_in = (p,l) -> callback_global(p,l,data,M,unique_t)
  
      # lb,ub = options
      bounds = options
  
      loss_in(pinit)
  
      # sr = [(0.0,1.0),(0.0,100.0)]
      # res_bbo = bboptimize(loss_in; SearchRange = sr,
      # res_bbo = bboptimize(loss_in; SearchRange = options,
      bounds_bbo = map(x-> (x[1], x[2]) , bounds)
      res_bbo = bboptimize(loss_in; SearchRange = bounds_bbo,
                                  NumDimensions = length(pinit),
                                  Method = :adaptive_de_rand_1_bin_radiuslimited,
                                  # Method = :separable_nes,
                                  NThreads=Threads.nthreads()-1,
                                  MaxSteps = 10000)#,
  
      p₀_bbo = best_candidate(res_bbo)
      @show p₀_bbo

      # may uncomment the the following line to feed initial guess for optimization
      p₀_bbo = pinit
       
      lb,ub = boundaries(bounds)
      res = optimize(loss_in,lb,ub,p₀_bbo, Fminbox(BFGS()), autodiff = :forward,
      # res = optimize(loss_in,lb,ub,p₀_bbo, Fminbox(NelderMead()), 
                                  Optim.Options(show_trace=true,
                                                f_tol = 1e-8, f_reltol=1e-5))#,
                                                # outer_iterations = 10))
  
      rates = res.minimizer
      l = res.minimum
      # bootstrap_rates = bootstrapping(rates,data,model,dens,cons,t
      hes = ForwardDiff.hessian(loss_in,rates)
      #finite difference
      # hes = FiniteDiff.finite_difference_hessian(loss_in,rates)
      println("Model used == ", M.model)
      println("Loss (MSE) == ", l)
      println("Rates == ", rates)
      res = OptRes(rates, l, hes)
      return res
end

using LinearAlgebra
# hessian2σ(H) = sqrt.(abs.( diag(inv(H)) ))
hessian2σ(H) = sqrt.( diag(inv(H)) )


function hessian2σ(H,L,data,rates)
  m = size(data,1)
  n = size(rates,1)
  ν = m - n # degrees of freedom
  MSE = L/ν # error nomralized by degrees of freedom
  # @show MSE
  σ² = abs.( diag(inv(H)) )
  # σ² = diag(inv(H))
  # @show σ²
  σ²_norm = σ² * MSE
  σ = sqrt.(σ²_norm)
  # σ = sqrt.(σ²)

  # σ² = diag( inv(H)*MSE ) # gives same results
  # σ²_norm = σ²
  # σ = sqrt.(σ²_norm)
  return σ
end

# function bootstrapping(rates,data,model,dens,cons,t)
#   #create sample
#   data = n
#   B_size = 10
#   # rand_mat = rand(1:length(data),(length(data),B_size))
#
#   rand_mat = [rand(1:length(data),length(data)) for i in 1:B_size]
#
#   bootstrap_rates = zeros(size(rates))
#   # bootstrap_rates = []
#   for i in 1:B_size
#     resampled_data = data[rand_mat[i]]
#     resampled_t = t[rand_mat[i]]
#     sort_ind = sortperm(resampled_t)
#     resampled_data = resampled_data[sort_ind]
#     resampled_t = resampled_t[sort_ind]
#     unique_t = unique(resampled_t)
#
#     # loss_in = (p) -> loss_2(p,resampled_data,model,dens,cons,resampled_t,unique_t)
#     # callback_in = (p,l) -> callback(p,l,resampled_data,model,dens,cons,resampled_t,unique_t)
#     loss_in = (p) -> loss(p,data,model,dens,cons,t,unique(t))
#     callback_in = (p,l) -> callback(p,l,data,model,dens,cons,t,unique(t))
#
#     # @show loss_in(rates)
#     res = optimize(loss_in,-9, 2, rates,Fminbox(NelderMead()))
#     # @time res = DiffEqFlux.sciml_train(loss_in,rates,BFGS(),
#     #                                 lower_bounds = -9,
#     #                                 upper_bounds = 2,
#     #                                 cb = callback_in,
#     #                                 maxiters=50,
#     #                                 f_tol = 1e-4)
#
#     # @show res.minimizer
#     bootstrap_rates= hcat(bootstrap_rates,res.minimizer)
#     # push!(bootstrap_rates,res.minimizer)
#     # Σ_sol = run_model_2(model,p,resampled_data,dens,cons,unique_t)
#   end
#
#   return bootstrap_rates[:,2:end]
# end
#

# function scan_param_space(p,data,model,dens,cons,t,unique_t,lb,ub)
#       #create rand mesh
#       n = 5000
#       Δbounds = ub - lb
#       # ps = rand(n,size(p,1)) .* 10 .- 8
#       ps = rand(n,size(p,1)) .* Δbounds  .+ lb
#       # ps = rand(n,size(p,1)) .* 100
#       errors = zeros(n)
#       for i in 1:n
#           # error, a = loss(ps[i,:],data,model,dens,cons,t)
#           error = loss(ps[i,:],data,model,dens,cons,t,unique_t)
#           errors[i] = error
#           # @show error
#           if mod(i,100) == 0
#             @show i, error
#           end
#       end
#       min_ind = argmin(errors)
#       min_error = minimum(errors)
#       @show findmin(errors)
#       p_best = ps[min_ind,:]
#       @show p_best
#       if size(p,1) == 2
#           fig = scatter(ps[:,1],ps[:,2],log10.(errors),
#           zcolor = log10.(errors),camera=[0,90])
#           scatter!([p_best[1]],[p_best[2]],[log10(min_error)],markercolor=:green,markersize=5)
#           display(fig)
#       end
#       return p_best
# end
