L1 = 3.700760241587129
L2 = 1.914718101631063
L3 = 2.4008921528083533
L4 = 36.00905301970478

Ls = [L1,L2,L3,L4]
nps = [1,3,3,3]
nps = [1,3,2,2]

AIC, w, rel_l = AIC_analysis(Ls,nps, 97)


function AIC_analysis(RSS,k,n)
    

    AICc = AIC_corrected(RSS, k, n);

    delta_i = AICc .- minimum(AICc)
    likelihoods = exp.(-delta_i ./ 2);
    sum_l = sum(likelihoods);
    w = likelihoods ./ sum_l; #this is w_i

    rel_l_mat = (w' * 1 ./ w)';
    # rel_l_mat = triu(rel_l_mat);
    return AICc, w, rel_l_mat      
end

function AIC_corrected(RSS, k, n)
    AIC = n .* log.(RSS ./ n) + 2 .* k; # AIC, where MLE of variance is RSS/n
    correction = 2 .* k .* (k .+ 1) ./ (n .- k .-1); # correction for small sample size
    AICc = AIC .+ correction;
end