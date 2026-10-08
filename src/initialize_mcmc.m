function vv = initialize_mcmc(N, T, M, ref_duration)
%INITIALIZE_MCMC Initialize binary raster with fixed spike count and refractory constraint.
%
% The function first draws M random occupied sites and then repairs
% refractory violations while preserving exactly M spikes.

    if M > N*T
        error('initialize_mcmc: M=%d exceeds total number of slots N*T=%d.', M, N*T);
    end

    vv = zeros(N, T);
    vv(randperm(N*T, M)) = 1;

    [vv_ref, ~, ~] = ref_period(vv, ref_duration);
    forbidden = find(vv & vv_ref);

    while ~isempty(forbidden)

        % Remove spikes that violate the refractory constraint.
        vv(forbidden) = 0;
        nf = numel(forbidden);

        % Recompute truly available slots after removing forbidden spikes.
        [~, ~, available_slots] = ref_period(vv, ref_duration);

        na = numel(available_slots);
        if na < nf
            error('initialize_mcmc: not enough available slots to preserve fixed spike count.');
        end

        % Reinsert exactly the number of removed spikes.
        vv(available_slots(randperm(na, nf))) = 1;

        [vv_ref, ~, ~] = ref_period(vv, ref_duration);
        forbidden = find(vv & vv_ref);
    end

    validate_state(vv, M, ref_duration);

end
