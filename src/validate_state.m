function validate_state(vv, M_expected, ref_duration)
%VALIDATE_STATE Check fixed spike count and refractory admissibility.
%
% vv is an N x T binary raster.
% M_expected is the required total number of spikes.
% ref_duration controls the absolute refractory constraint.

    if nnz(vv) ~= M_expected
        error('validate_state: spike count changed: nnz=%d, expected=%d.', ...
            nnz(vv), M_expected);
    end

    [vv_ref1, ~, ~] = ref_period(vv, ref_duration);

    if any(vv(:) & vv_ref1(:))
        error('validate_state: refractory constraint violated.');
    end

end
