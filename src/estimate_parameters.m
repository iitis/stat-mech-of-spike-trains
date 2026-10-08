function [a, b] = estimate_parameters(m1, sigma)
    m2 = sigma^2 + m1^2;

    x = sqrt(2*pi) * sigma * exp((m1^2) / (2*sigma^2));

    for k = 1:1000
        x_temp = x;
        x = quad(@(t) fun(t, x, m1, m2, sigma), 0, 1e10);
        if abs(x - x_temp) < 1e-10
            break
        end
    end

    a = (m2 - m1.*x) ./ ((sigma.^2) .* x);
    b = (x - m1) ./ (2 .* (sigma.^2) .* x);
end
