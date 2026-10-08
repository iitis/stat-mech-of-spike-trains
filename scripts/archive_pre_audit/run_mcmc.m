function run_mcmc(varargin)
% run_mcmc PREFIX N P_CONN DELMEAN DELSD T N_SPIKES REF_DURATION BETA N_ITER PHI SNAPSHOT_EVERY [SEED]
%
% Example:
%   run_mcmc('G',100,0.2,5,1,1000,500,2,10,20000,'tanh',2000,123)

    % ---- defaults ----
    prefix = 'G';
    N = 100;
    p_conn = 0.2;
    delmean = 5;
    delsd = 1;

    T = 100;
    n_spikes = 200;
    ref_duration = 3;
    beta = 1.0;
    n_iter = 20000;
    phi_name = 'tanh';
    snapshot_every = 5000;
    seed = []; % brak = losowo

    % ---- parsing of positional arguments ----
    if nargin >= 1, prefix = char(varargin{1}); end
    if nargin >= 2, N = str2double_if_needed(varargin{2}); end
    if nargin >= 3, p_conn = str2double_if_needed(varargin{3}); end
    if nargin >= 4, delmean = str2double_if_needed(varargin{4}); end
    if nargin >= 5, delsd = str2double_if_needed(varargin{5}); end
    if nargin >= 6, T = str2double_if_needed(varargin{6}); end
    if nargin >= 7, n_spikes = str2double_if_needed(varargin{7}); end
    if nargin >= 8, ref_duration = str2double_if_needed(varargin{8}); end
    if nargin >= 9, beta = str2double_if_needed(varargin{9}); end
    if nargin >=10, n_iter = str2double_if_needed(varargin{10}); end
    if nargin >=11, phi_name = char(varargin{11}); end
    if nargin >=12, snapshot_every = str2double_if_needed(varargin{12}); end
    if nargin >=13, seed = str2double_if_needed(varargin{13}); end

    if ~isempty(seed) && isfinite(seed)
        rng(seed, 'twister');
    end

    fprintf('PREFIX=%s N=%d p=%g delmean=%g delsd=%g T=%d spikes=%d ref=%d beta=%g n_iter=%d phi=%s snap=%d seed=%s\n', ...
        prefix, N, p_conn, delmean, delsd, T, n_spikes, ref_duration, beta, n_iter, phi_name, snapshot_every, mat2str(seed));

    % ---- choose nonlinearity ----
    phi = make_phi(phi_name);

    % ---- generate network ----
    G = initialize_network(N, p_conn, delmean, delsd);

    % ---- init microstate ----
    microstate = initialize_mcmc(G.N, T, n_spikes, ref_duration);

    % ---- operator ----
    B = transfer_operator(G.Adj, G.delays, T);

    % ---- MCMC ----
    [microstate, E_mem] = mcmc_core(microstate, G, B, beta, n_iter, ref_duration, phi, prefix, snapshot_every);

    % ---- write energy ----
    writematrix(E_mem, sprintf('%s_E_mem.txt', prefix));

    % ---- write final state ----
    write_microstate(microstate, sprintf('%s_microstate_iter_%d.txt', prefix, n_iter));

    % ---- save network params ----
    save(sprintf('%s_graph.mat', prefix), 'G');

    fprintf('Done.\n');
end

function x = str2double_if_needed(v)
    if ischar(v) || isstring(v)
        x = str2double(v);
    else
        x = v;
    end
end
