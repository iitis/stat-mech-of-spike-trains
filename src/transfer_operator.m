function B = transfer_operator(A, D, K)
%TRANSFER_OPERATOR Build temporal Fourier transfer operator for delayed propagation.
%
% A(i,j) = 1 denotes a directed edge i -> j.
% D(i,j) is the corresponding delay in time bins.
%
% MATLAB's fft convention is
%   n_hat(q) = sum_l n(l) exp(-2*pi*i*q*l/K).
%
% With this convention, a causal delay s(t - Delta) corresponds to
% multiplication by exp(-2*pi*i*q*Delta/K). We therefore use the negative
% phase explicitly and a non-conjugating transpose .'. This avoids relying
% on MATLAB's conjugating transpose operator ' to flip the sign implicitly.

    tt = 0:K-1;
    h = exp(-tt./5) - 0.5 .* exp(-(K-1-tt)./5);
    kern = fft(h);

    B = cell(1, K);

    for k = 1:K
        q = k - 1;

        % Causal delay phase.
        phase = exp(-2 .* pi .* 1i .* q ./ K .* D);

        % A(i,j) encodes i -> j, while the operator acts as
        % u_hat_j(q) = sum_i B_ji(q) n_hat_i(q).
        % The non-conjugating transpose .' performs the index swap only.
        B{k} = sparse(kern(k) .* (A .* phase).');
    end

end
