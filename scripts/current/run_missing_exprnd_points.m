repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(repo_root);
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'scripts', 'current')));

out_root = fullfile(repo_root, 'results', 'raw', 'missing_exprnd_fix');
if ~exist(out_root, 'dir')
    mkdir(out_root);
end

% Columns:
% label, beta, sdel, seed
jobs = {
    'b0_beta0p2_sdel4p5', 0.2, 4.5, 100201
    'b0_beta1p8_sdel4p5', 1.8, 4.5, 101801
    'b0_beta2p0_sdel4p5', 2.0, 4.5, 102001

    'b2_beta2p4_sdel4p5', 2.4, 4.5, 102401
    'b2_beta4p0_sdel4p5', 4.0, 4.5, 104001

    'b4_beta4p4_sdel4p5', 4.4, 4.5, 104401
    'b4_beta6p0_sdel4p5', 6.0, 4.5, 106001
};

for k = 1:size(jobs,1)
    label = jobs{k,1};
    beta  = jobs{k,2};
    sdel  = jobs{k,3};
    seed  = jobs{k,4};

    prefix = fullfile(out_root, label);

    fprintf(1, '\n=== MISSING POINT %d/%d: %s beta=%g sdel=%g seed=%d ===\n', ...
        k, size(jobs,1), label, beta, sdel, seed);

    run_phase_diagram( ...
        prefix, ...
        100, ...      % N
        0.3, ...      % p_conn
        5, ...        % delmean
        sdel, sdel, 0.15, ...   % one-point sdel grid
        500, ...      % T
        500, ...      % n_spikes
        3, ...        % ref_duration
        beta, beta, 0.2, ...    % one-point beta grid
        100000, ...   % n_iter
        10000, ...    % n_snap
        'heaviside3', ...
        seed ...
    );
end

fprintf(1, '\nAll missing-point reruns finished.\n');
