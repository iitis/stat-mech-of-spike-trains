function test_transfer_delay_sign()
%TEST_TRANSFER_DELAY_SIGN Check that transfer_operator implements causal delay.
%
% For an edge 1 -> 2 with delay Delta, a spike in neuron 1 at bin t0
% should produce the dominant response in neuron 2 around bin t0 + Delta
% modulo K.

    repo_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(repo_root, 'src')));

    N = 2;
    K = 16;
    delay = 3;

    A = zeros(N,N);
    D = zeros(N,N);

    A(1,2) = 1;      % edge 1 -> 2
    D(1,2) = delay;

    B = transfer_operator(A, D, K);

    v = zeros(N,K);
    t0 = 2;          % MATLAB 1-based time-bin index
    v(1,t0) = 1;

    [~, u] = apply_operator_nonlinear(v, B, @(x) x);

    [~, peak_bin] = max(u(2,:));
    expected_bin = mod((t0 - 1) + delay, K) + 1;

    fprintf('Presynaptic spike bin: %d\n', t0);
    fprintf('Delay: %d\n', delay);
    fprintf('Expected postsynaptic peak bin: %d\n', expected_bin);
    fprintf('Observed postsynaptic peak bin: %d\n', peak_bin);

    if peak_bin ~= expected_bin
        error('test_transfer_delay_sign: wrong delay direction. Expected bin %d, got bin %d.', ...
            expected_bin, peak_bin);
    end

    fprintf('test_transfer_delay_sign OK\n');

end
