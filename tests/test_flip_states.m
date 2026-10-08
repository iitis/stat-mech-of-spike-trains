function test_flip_states()
%TEST_FLIP_STATES Check that relocation moves preserve admissibility.

    repo_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(repo_root, 'src')));

    rng(1, 'twister');

    N = 20;
    T = 50;
    M = 120;
    ref_duration = 3;

    microstate = initialize_mcmc(N, T, M, ref_duration);
    validate_state(microstate, M, ref_duration);

    n_steps = 10000;

    for k = 1:n_steps
        microstate = flip_states(microstate, ref_duration);
        validate_state(microstate, M, ref_duration);
    end

    fprintf('test_flip_states OK\n');

end
