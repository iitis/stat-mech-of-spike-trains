repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'tests')));

fprintf('Running audit tests...\n');

test_transfer_delay_sign;
test_transfer_operator_old_new_equivalence;
test_initialize_mcmc;
test_flip_states;
test_walk_states;
test_move_invariants;
test_mcmc_core_xt;
test_mcmc_core_xtsnap;

fprintf('All enabled tests passed.\n');
