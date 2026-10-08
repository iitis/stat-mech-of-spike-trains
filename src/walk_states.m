function microstate_new = walk_states(microstate, ref_duration)
%WALK_STATES Symmetric single-attempt local temporal walk proposal.
%
% One existing spike is chosen uniformly and proposed to move by one time
% bin, either forward or backward with equal probability. If the target site
% is occupied or the resulting raster violates the refractory constraint,
% the original state is returned.
%
% This move preserves the total spike count and gives a symmetric proposal
% on the admissible state space, so it can be used with the standard
% Metropolis acceptance rule.

    microstate_new = microstate;

    occupied = find(microstate);
    if isempty(occupied)
        return
    end

    idx = occupied(randi(numel(occupied)));          % 1/M, symmetric
    [i, j] = ind2sub(size(microstate), idx);

    if rand < 0.5                                   % 1/2, symmetric
        j_new = modd(j + 1, size(microstate, 2));
    else
        j_new = modd(j - 1, size(microstate, 2));
    end

    if microstate(i, j_new)
        return
    end

    candidate = microstate;
    candidate(i, j) = 0;
    candidate(i, j_new) = 1;

    [vv_ref1, ~, ~] = ref_period(candidate, ref_duration);

    if ~any(candidate(:) & vv_ref1(:))
        microstate_new = candidate;
    end

end
