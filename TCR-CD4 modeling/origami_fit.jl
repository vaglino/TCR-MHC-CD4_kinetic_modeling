using Revise
using Optim
# change main directory accordingly ro where software is saved on own machine
main_dir = "C:/Users/stravaglino3/Documents/TCR-MHC-CD4 modeling for github/"

src_dir = main_dir*"src/"
includet(src_dir*"helper_functions.jl")

##Constants
mTCR = 25.0                         #ballpark coating density (um**-2) for TCR on the DNA Origami from Muaz's email
mTCR = mCD4 = 30.6 #±7.54                         #ballpark coating density (um**-2) for TCR on the DNA Origami from Muaz's email
                                     #For Origami there should be the same amount of CD4 as TCR
mpMHC = 550.0                         #ballpark coating density for MHC on the Origami target cell from Muaz's email

Pa_TCR_∞ = 0.1708                                 #Average Pa of the TCR Only Origami
n_TCR = Pa2n(Pa_TCR_∞)                     #average number of bonds
AcKₐTCR = n_TCR / (mTCR * mpMHC)              #Calculated Trimolecular 2D Affinity

Pa_CD4_∞ = 0.0142                                 #Average Pa of the CD4 Only Origami
n_CD4 = Pa2n(Pa_CD4_∞)                      #average number of bonds
AcKₐCD4 = n_CD4 / ((mCD4) * mpMHC)              #Calculated Trimolecular 2D Affinity

Pa_tot_∞ = 0.4159                                 #Average Total Pa of the Origami at 6nm
n_tot = Pa2n(Pa_tot_∞)                    #average number of total bonds at x=6nm
n_tri = n_tot - (n_CD4 + n_TCR)         #Average trimolecular number of bonds by subtracting the bimolecular from the total
AcKₐTri = n_tri / ((mCD4 * mTCR) * mpMHC)         #Calculated Trimolecular 2D Affinity

Ac = 3                                                        #Estimated contact area (um**2) for RBCs in Muaz's dissertation as an upper limit
max_separation = (2 * sqrt(Ac / π)) * 1000   #why x1000???    #This is the theoretical maximum separation distance (nm) in the contact area which is used as an x-value for the TCR-only condition

x_dis = [6., 13., 20., 100., max_separation]            #separation distances for the Origami
# x0 = 6.0                                                                #Offset for shortest separation distance

f₆ = 0.49                           #TCR bond fraction from the 6nm stiffness histogram data
f₁₃ = 0.62                          #TCR bond fraction from the 13nm stiffness histogram data
##

main_dir = "C:/Users/stravaglino3/Documents/TCR-MHC-CD4 modeling for github/"

data = CSV.read(main_dir*"TCR-CD4 modeling/Origami_data/Origami Normalized Affinity full data.csv",DataFrame)     #path to find the .csv data file

data_excluding_CD4only = data[!,:1:4] # get TCR-CD4-MHC data at different spacings

Y_data = delete_missing(data_excluding_CD4only)

Y = reduce(vcat,Y_data) 
σ_data = std.(Y_data)

## ORGANIZE DATA
X = []
for i in 1:length(Y_data)
    X = vcat( X, fill(x_dis[i], size(Y_data[i])) )
end
X

σ = []
for i in 1:length(Y_data)
    σ = vcat( σ, fill(σ_data[i], size(Y_data[i])) )
end
σ
##

function origami_model(p,X)
    KₐCD4,x₀,l = p
    n_tot_norm = AcKₐTCR * ( 1 + mCD4*KₐCD4*exp( -(X-x₀)/l ) )
end

# CHI SQUARE LOSS 
function χ²(Y,yₘ,σ)
    χ² = sum( ((Y .- yₘ)./σ ).^2 )
    # χ² = sum( ((Y .- yₘ) ).^2 )

end

# λ₆(Y,X,yₘ,σ,l) = sum( 2*AcKₐTCR .* (Y .- yₘ) ./σ.^2 .* exp.((6 .- X)./l) .* (13 .- X)./7 )
# λ₁₃(Y,X,yₘ,σ,l) = sum( 2*AcKₐTCR .* (Y .- yₘ) ./σ.^2 .* exp.((13 .- X)./l) .* (X .- 6)./7 )

function g(p,f,sep)
    KₐCD4,x₀,l = p
    g = 1 - 1/f + mCD4*KₐCD4*exp(-(sep-x₀)/l)
end

## loss function wrapper
function L(params,Y,X,σ)

    p = params

    yₘ = map(X -> origami_model(p,X), X)

    χ2 = χ²(Y,yₘ,σ)
    L = χ2 
    # @show L
    return L
end
    

# #Minimize the Cost Function
# x_0 = 6                                     #Initial Guesses for all the variables used to minimize the cost function
# l = 5
AcKₐTCR = mean(skipmissing(data[!,5])) 

Kₐ_CD4 = AcKₐCD4 / Ac       #Intial guess for the affinity of CD4 binding to TCR-MHC is set equal to CD4 binding to pMHC alone affinity
# lambda6 = 1
# lambda13 = 1 
              # Kₐ_CD4,  x₀, l, λ6, λ13  
# initial_guess = [Kₐ_CD4, 6., 5., 0., 0.]
initial_guess = [Kₐ_CD4, 6., 5.] 

# bounds = [(0., Inf), (0, Inf), (0, Inf), (-Inf,Inf), (-Inf,Inf)]#, (0, Inf), (0, Inf)]
bounds = [(0., 100), (0, 100), (0, 100)]#, (0, Inf), (0, Inf)]
# bounds = [(0., 100), (0, 100), (0, 100), (0, 10), (0, 10)]


# OPTIMIZATION ROUTINE
lb,ub = boundaries(bounds)
L_optimization = (p) -> L(p,Y,X,σ)
@show L(initial_guess,Y,X,σ)

# using BlackBoxOptim
# res_bbo = bboptimize(L_optimization; SearchRange = bounds,
#                                   NumDimensions = length(initial_guess),
#                                   Method = :adaptive_de_rand_1_bin_radiuslimited,
#                                   # Method = :separable_nes,
#                                   NThreads=Threads.nthreads()-1,
#                                   MaxSteps = 100000)#,

# initial_guess = best_candidate(res_bbo)
od = OnceDifferentiable(L_optimization, initial_guess; autodiff = :forward);
res = optimize(od,lb,ub,initial_guess,Fminbox(BFGS()) ,
# res = optimize( od,initial_guess,BFGS() )#,
# res = optimize(L_optimization,lb,ub,p₀_bbo,Fminbox(NelderMead()),#
                                Optim.Options(show_trace=true,
                                f_tol = 1e-5, f_reltol=1e-5))


L(res.minimizer,Y,X,σ)

# using FiniteDiff
# # wres = sp.minimize(chi2, initial_guess, bounds = bounds )#, method = "trust-constr")         #This is the minimization function scipy.optimize.minimize() with the function to minimize, initial guesses, boundaries on the variables, and a commented out option to change the algorithm used (method)
# @show res.minimizer 
# @show res.minimum 
# # hes = ForwardDiff.hessian(L_optimization,initial_guess)
# hes = ForwardDiff.hessian(L_optimization,res.minimizer)

# # hes = FiniteDiff.finite_difference_hessian(L_optimization,res.minimizer)
# σ_params = hessian2σ(hes) 
# hessian2σ(H) = sqrt.(diag(inv(H)))
# inv(hes)
# using Measurements

# params = res.minimizer .± σ_params
# params = res.minimizer
# ##
# xrange = 5:100
# p1 = scatter(X,Y)
# ymodel = map(x -> origami_model(params,x), xrange)
# # ymodel = map(x -> origami_model(initial_guess,x), xrange)

# plot!(xrange,ymodel)
# plotly()
# p2 = scatter(log10.(X),Y)
# plot!(log10.(xrange),ymodel)

# plot(p1,p2)
# ##

# CSV.write(res_dir*"origami_fit.csv", measurements2Df(xrange,ymodel;labels=["fit","std"]))    

# BOOTSTRAPPING FUNCTION 
function bootstrap_fit(p,X,Y; boot_n=1000)
    n_data = length(Y)
    p_boot = []
    for n in 1:boot_n
        sampling = rand(1:n_data, n_data)
        Xₙ = X[sampling]
        Yₙ = Y[sampling]
        σₙ = σ[sampling]
        L_optimization = (p) -> L(p,Yₙ,Xₙ,σₙ)
        od = OnceDifferentiable(L_optimization, p; autodiff = :forward);
        res = optimize(od,lb,ub,p,Fminbox(BFGS()),
                                    Optim.Options(show_trace=false,
                                    f_tol = 1e-8, f_reltol=1e-8))
        push!(p_boot,res.minimizer)
    end
    mean(reduce(vcat,p_boot),dims=1)
    return p_boot
end

p_boot = bootstrap_fit(res.minimizer,X,Y;boot_n=1000)
params_boot = mean(p_boot) .± std(p_boot)


## USE AcKaTCR with errors for propagating during plotting
AcKₐTCR = mean(skipmissing(data[!,5])) .± std(skipmissing(data[!,5]))

## PLOT BOTH IN LINEAR AND LOG SPACE
plotly()
xrange = 5:120
ymodel = map(x -> origami_model(params_boot,x), xrange)

p1 = scatter(X,Y)
plot!(xrange,ymodel)

p2 = scatter(log10.(X),Y)
plot!(log10.(xrange),ymodel)

plot(p1,p2)
##

# uncomment to save results
# CSV.write(res_dir*"origami_fit.csv", measurements2Df(xrange,ymodel;labels=["fit","std"])) 
