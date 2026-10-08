function [vv_ref1, vv_ref2, available_slots] = ref_period(vv, ref_duration)

    vv_ref1 = ref_period_convolve(vv, ref_duration);
    vv_ref2 = vv | vv_ref1 | fliplr(ref_period_convolve(fliplr(vv), ref_duration));
    available_slots = find(~vv_ref2);

end
