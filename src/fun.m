function y = fun(t, x, m1, m2, sigma)
    a = (m2 - m1.*x) ./ ((sigma.^2) .* x);
    b = (x - m1) ./ (2 .* (sigma.^2) .* x);
    y = exp(-a.*t - b.*(t.^2));
end
