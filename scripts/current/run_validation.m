% Validate the sampler and historical reconstruction without production runs.
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(repo_root);
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'tests')));
fprintf('MATLAB: %s\nPlatform: %s\nRepository: %s\n', version, computer, repo_root);

fprintf('\n=== 1/6: sampler tests ===\n');
run(fullfile(repo_root, 'tests', 'run_all_tests.m'));
fprintf('\n=== 2/6: smoke simulation ===\n');
run(fullfile(repo_root, 'scripts', 'current', 'run_phase_smoke_test.m'));
fprintf('\n=== 3/6: reconstruct historical main map ===\n');
run(fullfile(repo_root, 'scripts', 'current', 'merge_missing_exprnd_points.m'));
fprintf('\n=== 4/6: verify reconstruction and exports ===\n');
test_main_map_reconstruction;
fprintf('\n=== 5/6: repeat reconstruction and verify again ===\n');
run(fullfile(repo_root, 'scripts', 'current', 'merge_missing_exprnd_points.m'));
test_main_map_reconstruction;
fprintf('\n=== 6/6: regenerate entropy/IPR figures ===\n');
set(groot, 'defaultFigureVisible', 'off');
run(fullfile(repo_root, 'scripts', 'current', 'plot_entropy_maps_matlab.m'));
fprintf('\nVALIDATION PASSED: sampler tests, smoke run, reconstruction, exports, repeated reconstruction, figures.\n');
