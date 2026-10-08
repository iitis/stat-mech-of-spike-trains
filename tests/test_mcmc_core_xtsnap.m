function test_mcmc_core_xtsnap()
%TEST_MCMC_CORE_XTSNAP Check snapshot MCMC core and diagnostics.

    repo_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(repo_root, 'src')));

    rng(5, 'twister');

    N = 12;
    T = 40;
    M = 50;
    ref_duration = 3;
    beta = 1.0;
    n_iter = 200;
    n_snap = 50;

    p_conn = 0.2;
    delmean = 5;
    delsd = 1;

    G = initialize_network(N, p_conn, delmean, delsd);
    B = transfer_operator(G.Adj, G.delays, T);

    phi = @(x) tanh(x);

    microstate0 = initialize_mcmc(N, T, M, ref_duration);

    [x, Efin, snaps, diag] = mcmc_core_xtsnap(microstate0, G, B, beta, n_iter, ref_duration, phi, n_snap, M);

    if numel(x) ~= T
        error('test_mcmc_core_xtsnap: x has wrong length.');
    end

    if ~isfinite(Efin)
        error('test_mcmc_core_xtsnap: Efin is not finite.');
    end

    expected_snaps = floor(n_iter / n_snap);

    if size(snaps, 1) ~= expected_snaps
        error('test_mcmc_core_xtsnap: wrong number of snapshots. Expected %d, got %d.', ...
            expected_snaps, size(snaps, 1));
    end

    if ~diag.spike_count_ok
        error('test_mcmc_core_xtsnap: spike count diagnostic failed.');
    end

    if ~diag.refractory_ok
        error('test_mcmc_core_xtsnap: refractory diagnostic failed.');
    end

    fprintf('test_mcmc_core_xtsnap OK, accept_rate = %.3f, snaps = %d\n', ...
        diag.accept_rate, diag.snap_count);

end
