function sample = rounded_randn(mean0, std0, n)
    mean0 = mean0 - 0.5;
    sample = ceil(std0 .* randn(n,1) + mean0);

    std_sample = std(sample);
    std_corrected = std0;

    step = 0.01; temp1 = 0; temp2 = 0; k = 1;

    while abs(std0 - std_sample) > 0.001 && k < 10000
        if std_sample < std0
            std_corrected = std_corrected + step;
            if temp1 == -1
                step = 0.1 * step;
            end
            temp1 = 1;
            temp2 = temp2 + 1;
            if temp2 == 10
                temp2 = 0;
                step = 10 * step;
            end
        else
            std_corrected = std_corrected - step;
            if temp1 == 1
                step = 0.1 * step;
                temp2 = 0;
            end
            temp1 = -1;
        end

        sample = ceil(std_corrected .* randn(n,1) + mean0);
        std_sample = std(sample);
        k = k + 1;
    end

    % disp([k, std_corrected, std_sample, mean(sample)]);
end
