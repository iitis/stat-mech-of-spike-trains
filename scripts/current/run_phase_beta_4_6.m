repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(repo_root);
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'scripts', 'current')));

prefix = fullfile(repo_root, 'results', 'raw', 'beta_0_6_sparse_1M', ...
    'Gphase_sparse_sdel0to4p5_beta4to6_h3_1M');

run_phase_diagram( ...
    prefix, ...
    100, ...      % N
    0.3, ...      % p_conn
    5, ...        % delmean
    0.0, 4.5, 0.9, ...      % sdel grid: 6 points
    500, ...      % T
    500, ...      % n_spikes
    3, ...        % ref_duration
    4.5, 6.0, 0.5, ...      % beta grid: 4 points
    1000000, ...  % n_iter
    100000, ...   % n_snap
    'heaviside3', ...
    1 ...
);