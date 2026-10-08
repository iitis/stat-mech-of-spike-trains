function phi = make_phi(name)
    name = lower(string(name));

    switch name
        case "tanh"
            phi = @(x) tanh(x);
        case "relu"
            phi = @(x) max(0, x);
        case "heaviside6"
            phi = @(x) single(x >= 6);
        case "heaviside12"
            phi = @(x) single(x >= 12);
        case "heaviside3"
            phi = @(x) single(x >= 3);    
        case "sigmoid"
            phi = @(x) 1 ./ (1 + exp(-x));
        case "identity"
            phi = @(x) x;
        otherwise
            error("Unknown phi: %s (use tanh|relu|heaviside|sigmoid|identity)", name);
    end
end
