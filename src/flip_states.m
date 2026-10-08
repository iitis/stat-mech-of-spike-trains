function microstate_new = flip_states(microstate, ref_duration)
%FLIP_STATES Symmetric spike relocation proposal.
%
% One occupied site and one empty site are selected uniformly. The spike is
% relocated only if the resulting raster satisfies the refractory constraint.
% Otherwise the proposal returns the original state.
%
% This move preserves the total spike count and gives a symmetric proposal
% on the admissible state space, so it can be used with the standard
% Metropolis acceptance rule.

    microstate_new = microstate;

    occupied_slots = find(microstate);
    empty_slots = find(~microstate);

    if isempty(occupied_slots) || isempty(empty_slots)
        return
    end

    old_slot = occupied_slots(randi(numel(occupied_slots)));
    new_slot = empty_slots(randi(numel(empty_slots)));

    candidate = microstate;
    candidate(old_slot) = 0;
    candidate(new_slot) = 1;

    [vv_ref1, ~, ~] = ref_period(candidate, ref_duration);

    if ~any(candidate(:) & vv_ref1(:))
        microstate_new = candidate;
    end

end
