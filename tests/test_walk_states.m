function test_walk_states()
%TEST_WALK_STATES Check that local walk proposals preserve admissibility.

    repo_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(repo_root, 'src')));

    rng(2, 'twister');

    N = 20;
    T = 50;
    M = 120;
    ref_duration = 3;

    microstate = initialize_mcmc(N, T, M, ref_duration);
    validate_state(microstate, M, ref_duration);

    n_steps = 10000;
    n_changed = 0;

    for k = 1:n_steps
        old_state = microstate;
        microstate = walk_states(microstate, ref_duration);
        validate_state(microstate, M, ref_duration);

        if any(microstate(:) ~= old_state(:))
            n_changed = n_changed + 1;
        end
    end

    fprintf('test_walk_states OK, changed proposals: %d/%d\n', n_changed, n_steps);

end
