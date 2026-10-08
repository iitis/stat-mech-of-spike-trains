function test_initialize_mcmc()
%TEST_INITIALIZE_MCMC Check that initialization preserves fixed spike count
% and satisfies the refractory constraint.

    repo_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(repo_root, 'src')));

    N = 20;
    T = 50;
    M = 120;
    ref_duration = 3;

    for seed = 1:100
        rng(seed, 'twister');
        vv = initialize_mcmc(N, T, M, ref_duration);
        validate_state(vv, M, ref_duration);
    end

    fprintf('test_initialize_mcmc OK\n');

end
