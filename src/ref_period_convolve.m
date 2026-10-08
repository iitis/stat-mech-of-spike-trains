function vv_ref = ref_period_convolve(vv, ref_duration)

    vv_ref = zeros(size(vv,1), size(vv,2) + ref_duration);

    for k = 1:size(vv,1)
        vv_ref(k,:) = conv(vv(k,:), [0, ones(1, ref_duration)]);
    end

    % time-periodic boundary conditions
    T = size(vv,2);
    vv_ref(:,1:ref_duration) = vv_ref(:,1:ref_duration) | vv_ref(:,T+1:T+ref_duration);

    vv_ref = vv_ref(:,1:T);

end
