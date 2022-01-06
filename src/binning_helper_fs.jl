
delete_missing(x) = [collect(skipmissing(x[:,i])) for i in 1:ncol(x)]

function bin_dissociations(data, n_bins; bin_type = "force", edges=nothing)
    F = data[1]
    t = data[2]
    # make sure data is sorted by increasing F
    sort_i = sortperm(F)
    F = F[sort_i]
    t = t[sort_i]
    # even force bins
    if edges != nothing
        bins,bin_ids = partition(F,edges)
        binned_t = [similar(t,0) for i in bins]

    elseif bin_type == "force"
        f_min = minimum(F)
        f_max = maximum(F)

        edges = range(f_min,stop=f_max,length=n_bins+1)
        bins,bin_ids = partition(F,edges)
        binned_t = [similar(t,0) for i in bins]

    elseif bin_type == "equal number"
        edges_ind = Int64.(round.(range(1,stop=length(F),length=n_bins+1)))
        edges = F[edges_ind]
        # edges = [0,9,12,15,19,30]
        bins,bin_ids = partition(F,edges)
        binned_t = [similar(t,0) for i in bins]
    end
    # h = fit(Histogram, F)

    for (i,bin) in enumerate(bin_ids)
        push!(binned_t[bin],t[i])
    end
    
    return bins, binned_t, edges
end


function partition(data, intervals)
    ranges = intervals[1:end-1] .=> intervals[2:end]
    bins = [similar(data, 0) for _ in 1:length(ranges)]
    ids = similar(data,0)
    for x in data
        for (i, (a, b)) in pairs(ranges)
            if i != length(bins) && a <= x < b
                push!(bins[i], x)
                push!(ids,i)
                break

            elseif i == length(bins) && a <= x <= b
                push!(bins[i], x)
                push!(ids,i)
                break
        end
        end
    end
    return bins, Int64.(ids)
end

function resort(f,t)
    t_sorted = sort.(t)
    sort_id = sortperm.(t)
    f_sorted = [f[i][sort_id[i]] for i in 1:length(sort_id)]
    return f_sorted, t_sorted
end

survival_p(t) = collect((length(t) :-1:1) / length(t))

p2lnp(p) = map((x) -> log.(x), p)


function bootstrap_bins(data, n_bins, n_boot; bin_type="force",edges=nothing)
    # create bootstrap samples
    d = Normal()
    n_data = length(data[1])
    boot_samples = []
    for n in 1:n_boot 
        sampling = rand(1:n_data, n_data)
        σ = rand(Normal(0.,1), n_data)
        sampleₙ = [data[1][sampling] + σ ,data[2][sampling]]
        Fₙ, tₙ, edgesₙ = bin_dissociations(sampleₙ, n_bins; 
                        bin_type = bin_type,
                        edges=edges)
        bootₙ = Ft_average_and_error(Fₙ, tₙ)
        push!(boot_samples,bootₙ)
    end
    Fs, ts, edges = bin_dissociations(data, n_bins; 
                        bin_type = bin_type,
                        edges=edges)
    println("bin edges = ", edges)
    println("number of points per bin = ",map(F->length(F),Fs))

    f_means = []
    f_sems = []
    t_means = []
    t_sems = []
    for bootₙ in boot_samples
        push!(f_means,bootₙ.f_mean)
        push!(f_sems,bootₙ.f_sem)
        push!(t_means,bootₙ.t_mean)
        push!(t_sems,bootₙ.t_sem)
    end
    f_means = reduce(hcat,f_means) # reorganize from array of vectors to matrix
    f_sems = reduce(hcat,f_sems)
    t_means = reduce(hcat,t_means)
    t_sems = reduce(hcat,t_sems)

    f_mean = vec(nanmean(f_means,2))
    t_mean = vec(nanmean(t_means,2))
    f_sem = vec(nanmean(f_sems,2))
    t_sem = vec(nanmean(t_sems,2))
 
    # @show f_mean
    # f_mean = mean(f_means)
    # t_mean = mean(t_means)
    # f_sem = mean(f_sems)
    # t_sem = mean(t_sems)
    df = DataFrame(f_mean=f_mean,f_sem=f_sem,
                    t_mean=t_mean,t_sem=t_sem)
end


nanmean(x) = mean(filter(!isnan,x))
nanmean(x,dims) = mapslices(nanmean,x,dims=dims)

using Distributions
function combined_bootstrap_dataset(data, n_boot)
    # create bootstrap samples
    d = Normal()
    n_data = length(data[1])
    F = []
    t = []
    for n in 1:n_boot 
        
        sampling = rand(1:length(data[1]), length(data[1]))
        σ = rand(Normal(0., 0.3), n_data)
        # sampleₙ = [data[1][sampling],data[2][sampling]]
        Fₙ = data[1][sampling] + σ
        tₙ = data[2][sampling]
        # Fₙ, tₙ, edgesₙ = bin_dissociations(sampleₙ, n_bins; 
        #                 bin_type = bin_type,
        #                 edges=edges)
        # bootₙ = Ft_average_and_error(Fₙ, tₙ)
        # @show size(sampleₙ)
        # push!(boot_samples,sampleₙ )
        F = vcat(F,Fₙ)
        t = vcat(t,tₙ)
    end
    F = Float64.(F)
    t = Float64.(t)
    # combined_samples = reduce(hcat,boot_samples)
    return [F,t]
end
    

