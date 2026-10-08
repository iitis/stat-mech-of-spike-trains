repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'scripts', 'current')));

prefix = fullfile(repo_root, 'results', 'raw', 'smoke_test', 'Gphase_smoke');

run_phase_diagram( ...
    prefix, ...
    20, ...       % N
    0.2, ...      % p_conn
    5, ...        % delmean
    0.0, 0.4, 0.2, ...      % sdel grid
    50, ...       % T
    40, ...       % n_spikes
    3, ...        % ref_duration
    0.0, 0.4, 0.2, ...      % beta grid
    2000, ...     % n_iter
    1000, ...     % n_snap
    'tanh', ...
    123 ...
);
