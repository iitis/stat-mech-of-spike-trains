function s = me_sample(m1, sigma, n)
    [a, b] = estimate_parameters(m1, sigma);
    s = inv_cdf(rand(n,1), a, b);
end
