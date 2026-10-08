function test_move_invariants()
%TEST_MOVE_INVARIANTS Check invariants of walk and relocation moves.
%
% The local walk move should preserve per-neuron spike counts.
% The relocation move should be able to change per-neuron spike counts.

    repo_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(repo_root, 'src')));

    rng(3, 'twister');

    N = 20;
    T = 50;
    M = 120;
    ref_duration = 3;

    microstate = initialize_mcmc(N, T, M, ref_duration);
    validate_state(microstate, M, ref_duration);

    counts0 = sum(microstate, 2);

    % Walk moves should preserve per-neuron spike counts.
    x = microstate;
    for k = 1:1000
        x = walk_states(x, ref_duration);
        validate_state(x, M, ref_duration);

        counts = sum(x, 2);
        if any(counts ~= counts0)
            error('test_move_invariants: walk move changed per-neuron spike counts.');
        end
    end

    % Relocation moves should be capable of changing per-neuron counts.
    y = microstate;
    changed_counts = false;

    for k = 1:10000
        y = flip_states(y, ref_duration);
        validate_state(y, M, ref_duration);

        counts = sum(y, 2);
        if any(counts ~= counts0)
            changed_counts = true;
            break
        end
    end

    if ~changed_counts
        error('test_move_invariants: relocation did not change per-neuron counts in 10000 attempts.');
    end

    fprintf('test_move_invariants OK\n');

end
