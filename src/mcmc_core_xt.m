function [x, E_final, diag] = mcmc_core_xt(microstate, G, B, beta, n_iter, ref_duration, phi, M_expected)
%MCMC_CORE_XT Run Metropolis sampling and return final population activity.
%
% [x, E_final] = mcmc_core_xt(...)
% [x, E_final, diag] = mcmc_core_xt(..., M_expected)
%
% The function samples binary spike rasters using a mixture of symmetric
% single-attempt relocation and temporal-walk proposals.

    %#ok<NASGU> % G is kept for interface compatibility.

    T = size(microstate, 2);

    if nargin < 8 || isempty(M_expected)
        M_expected = nnz(microstate);
    end

    validate_state(microstate, M_expected, ref_duration);

    diag = struct();
    diag.n_iter = n_iter;
    diag.M_expected = M_expected;
    diag.initial_spike_count = nnz(microstate);
    diag.n_accept = 0;
    diag.n_reject = 0;
    diag.n_flip = 0;
    diag.n_walk = 0;

    % ---- initial energy ----
    E = -apply_operator_nonlinear(microstate, B, phi);

    for k = 1:n_iter

        % ---- propose ----
        if rand < 0.5
            microstate_new = flip_states(microstate, ref_duration);
            diag.n_flip = diag.n_flip + 1;
        else
            microstate_new = walk_states(microstate, ref_duration);
            diag.n_walk = diag.n_walk + 1;
        end

        % Audit safety: proposals should preserve the ensemble constraints.
        if mod(k, 1000) == 0
            validate_state(microstate_new, M_expected, ref_duration);
        end

        % ---- evaluate ----
        E_new = -apply_operator_nonlinear(microstate_new, B, phi);

        % ---- accept/reject ----
        acc = min(1, exp(-beta * (E_new - E)));

        if rand < acc
            microstate = microstate_new;
            E = E_new;
            diag.n_accept = diag.n_accept + 1;
        else
            diag.n_reject = diag.n_reject + 1;
        end
    end

    validate_state(microstate, M_expected, ref_duration);

    % ---- output order signal x(t) ----
    x = sum(double(microstate), 1);
    x = x(:).';

    if numel(x) ~= T
        error('mcmc_core_xt: x has wrong length (expected %d).', T);
    end

    E_final = E;

    diag.final_spike_count = nnz(microstate);
    diag.spike_count_ok = (diag.final_spike_count == M_expected);

    [vv_ref1, ~, ~] = ref_period(microstate, ref_duration);
    diag.refractory_ok = ~any(microstate(:) & vv_ref1(:));

    diag.accept_rate = diag.n_accept / max(1, n_iter);

    if ~diag.spike_count_ok
        error('mcmc_core_xt: final spike count changed: nnz=%d expected=%d.', ...
            diag.final_spike_count, M_expected);
    end

    if ~diag.refractory_ok
        error('mcmc_core_xt: final refractory constraint violated.');
    end

end
