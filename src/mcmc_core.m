function [microstate, E_mem] = mcmc_core(microstate, G, B, beta, n_iter, ref_duration, phi, prefix, snapshot_every)

    E_mem = NaN(n_iter,1);

    % ---- initial energy ----
    E = -apply_operator_nonlinear(microstate, B, phi);

    for k = 1:n_iter
        if rand < 0.5
            microstate_new = flip_states(microstate, ref_duration);
        else
            microstate_new = walk_states(microstate, ref_duration);
        end

        E_new = -apply_operator_nonlinear(microstate_new, B, phi);

        if rand < exp(-beta*(E_new - E))
            microstate = microstate_new;
            E = E_new;
        end

        E_mem(k) = E;

        % ---- log modulo 1000 ----
        if mod(k,1000)==0
            fprintf('iter %d E=%g\n', k, E);
        end

        % ---- snapshot ----
        if snapshot_every > 0 && mod(k, snapshot_every)==0
            write_microstate(microstate, sprintf('%s_microstate_iter_%d.txt', prefix, k));
        end
    end
end
