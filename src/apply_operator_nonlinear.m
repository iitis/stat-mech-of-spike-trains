function [E, v1] = apply_operator_nonlinear(v, B, phi)
    % v: N x K (0/1)
    % B: cell(1,K), B{k} sparse NxN
    % phi: function handle

    u1 = fft(double(v), [], 2);
    K = size(u1, 2);

    if length(B) ~= K
        error('apply_operator_nonlinear: length(B) ~= K');
    end

    for k = 1:K
        u1(:,k) = B{k} * u1(:,k);
    end

    v1 = real(ifft(u1, [], 2));   % N x K

    % apply nonlinearity
    E = sum(sum( double(v) .* phi(v1) ));
end
