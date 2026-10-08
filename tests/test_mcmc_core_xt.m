function test_mcmc_core_xt()
%TEST_MCMC_CORE_XT Check that MCMC core preserves ensemble constraints.

    repo_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(repo_root, 'src')));

    rng(4, 'twister');

    N = 12;
    T = 40;
    M = 50;
    ref_duration = 3;
    beta = 1.0;
    n_iter = 200;

    p_conn = 0.2;
    delmean = 5;
    delsd = 1;

    G = initialize_network(N, p_conn, delmean, delsd);
    B = transfer_operator(G.Adj, G.delays, T);

    phi = @(x) tanh(x);

    microstate0 = initialize_mcmc(N, T, M, ref_duration);

    [x, Efin, diag] = mcmc_core_xt(microstate0, G, B, beta, n_iter, ref_duration, phi, M);

    if numel(x) ~= T
        error('test_mcm_core_xt: x has wrong length.');
    end

    if ~isfinite(Efin)
        error('test_mcmc_core_xt: Efin is not finite.');
    end

    if ~diag.spike_count_ok
        error('test_mcmc_core_xt: spike count diagnostic failed.');
    end

    if ~diag.refractory_ok
        error('test_mcmc_core_xt: refractory diagnostic failed.');
    end

    fprintf('test_mcmc_core_xt OK, accept_rate = %.3f\n', diag.accept_rate);

end
