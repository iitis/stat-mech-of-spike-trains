function run_phase_diagram(varargin)
% run_phase_diagram PREFIX N P_CONN DELMEAN SDEL_MIN SDEL_MAX SDEL_STEP T N_SPIKES REF_DURATION ...
%                  BETA_MIN BETA_MAX BETA_STEP N_ITER N_SNAP PHI [SEED]
%
% Example:
% run_phase_diagram('G1',100,0.2,5, 0.0,2.0,0.2, 200,200,3, 0.2,2.0,0.2, 20000,100,'tanh',123)

    % ---- defaults ----
    prefix = 'Gphase';
    N = 100;
    p_conn = 0.2;
    delmean = 5;

    sdel_min = 0.0;
    sdel_max = 2.0;
    sdel_step = 0.2;

    T = 200;
    n_spikes = 200;
    ref_duration = 3;

    beta_min = 0.1;
    beta_max = 2.0;
    beta_step = 0.2;

    n_iter = 20000;
    n_snap = inf;      % <--- NEW: domyślnie brak snapshotów
    phi_name = 'tanh';
    seed = [];

    % ---- parse args ----
    if nargin >= 1,  prefix = char(varargin{1}); end
    if nargin >= 2,  N = str2double_if_needed(varargin{2}); end
    if nargin >= 3,  p_conn = str2double_if_needed(varargin{3}); end
    if nargin >= 4,  delmean = str2double_if_needed(varargin{4}); end
    if nargin >= 5,  sdel_min = str2double_if_needed(varargin{5}); end
    if nargin >= 6,  sdel_max = str2double_if_needed(varargin{6}); end
    if nargin >= 7,  sdel_step = str2double_if_needed(varargin{7}); end
    if nargin >= 8,  T = str2double_if_needed(varargin{8}); end
    if nargin >= 9,  n_spikes = str2double_if_needed(varargin{9}); end
    if nargin >= 10, ref_duration = str2double_if_needed(varargin{10}); end
    if nargin >= 11, beta_min = str2double_if_needed(varargin{11}); end
    if nargin >= 12, beta_max = str2double_if_needed(varargin{12}); end
    if nargin >= 13, beta_step = str2double_if_needed(varargin{13}); end
    if nargin >= 14, n_iter = str2double_if_needed(varargin{14}); end
    if nargin >= 15, n_snap = str2double_if_needed(varargin{15}); end   % <--- NEW
    if nargin >= 16, phi_name = char(varargin{16}); end                 % <--- shifted
    if nargin >= 17, seed = str2double_if_needed(varargin{17}); end     % <--- shifted

    % normalize n_snap
    if isempty(n_snap) || ~isfinite(n_snap) || n_snap <= 0
        n_snap = inf;
    else
        n_snap = floor(n_snap);
    end

    if ~isempty(seed) && isfinite(seed)
        rng(seed, 'twister');
    end

    % ---- grids ----
    betas = beta_min:beta_step:beta_max;
    sdels = sdel_min:sdel_step:sdel_max;

    nb = numel(betas);
    ns = numel(sdels);

    % ---- choose nonlinearity ----
    phi = make_phi(phi_name);

    % ---- output file ----
    out_file = sprintf('%s_phase.mat', prefix);
    if exist(out_file, 'file')
        delete(out_file);
    end

    % ---- debug metadata (plain char = less MATLAB v7.3 MCOS weirdness) ----
    code_version = 'run_phase_diagram DEBUG 2026-02-04 SNAPSHOT_V73';
    code_pwd = char(pwd);
    code_when = datestr(now);
    slurm_job_id = getenv('SLURM_JOB_ID');

    disp(code_version);
    disp(code_pwd);

    % ---- allocate in RAM (MCR-friendly) ----
    ok = false(nb, ns);
    exit_code = zeros(nb, ns, 'int32');
    nan_in_x = false(nb, ns);
    bad_len_x = false(nb, ns);

    x_grid = zeros(nb, ns, T, 'single');
    E_final = NaN(nb, ns);

    % ---- NEW: snapshot allocation ----
    do_snap = isfinite(n_snap);
    if do_snap
        max_snaps = floor(n_iter / n_snap);
    else
        max_snaps = 0;
    end

    % Snapshot: [nb ns max_snaps N T] uint8 (0/1)
    micro_snaps_grid = zeros(nb, ns, max_snaps, N, T, 'uint8');
    snap_count = zeros(nb, ns, 'int32'); % ile snapshotów realnie zapisano dla danego (ib,is)

    % ---- NEW: store initial network per grid point (compressed: edge list + delays) ----
    G_s_ind_grid = cell(nb, ns);        % uint32 linear indices of edges
    G_delay_vec_grid = cell(nb, ns);    % uint16/uint32 delays corresponding to s_ind
    G_E_grid = zeros(nb, ns, 'int32');  % number of edges
    G_D_grid = zeros(nb, ns, 'int32');  % max delay

    % ---- main loop ----
    fprintf(1, "SAVE_MODE = RAM_THEN_SAVE_V73\n");
    fprintf(1, "START %s pwd=%s out_file=%s nb=%d ns=%d T=%d SLURM_JOB_ID=%s\n", ...
            datestr(now), pwd, out_file, nb, ns, T, slurm_job_id);
    fprintf(1, "SNAPSHOT: do_snap=%d n_snap=%g max_snaps=%d\n", do_snap, n_snap, max_snaps);

    for ib = 1:nb
        beta = betas(ib);
        fprintf(1, "--- ENTER ib=%d/%d beta=%g ---\n", ib, nb, beta);

        for is = 1:ns
            delsd = sdels(is);

            fprintf(1, 'ib=%d/%d beta=%g  is=%d/%d sdel=%g ... ', ...
                    ib, nb, beta, is, ns, delsd);

            try
                G = initialize_network(N, p_conn, delmean, delsd);
                % ---- NEW: save initial network state for this (ib,is) ----
                % Store edges as linear indices + delays on those edges to avoid saving full NxN matrices.
                G_s_ind_grid{ib, is} = uint32(G.s_ind);

                if isempty(G.s_ind)
                    G_delay_vec_grid{ib, is} = uint16([]);
                else
                    dv = G.delays(G.s_ind);   % delays in the same order as s_ind (should be positive ints)
                    dv = dv(:);

                    % choose compact integer type
                    if max(dv) <= double(intmax('uint16'))
                        G_delay_vec_grid{ib, is} = uint16(dv);
                    else
                        G_delay_vec_grid{ib, is} = uint32(dv);
                    end
                end

                G_E_grid(ib, is) = int32(G.E);
                G_D_grid(ib, is) = int32(G.D);

                B = transfer_operator(G.Adj, G.delays, T);
                microstate0 = initialize_mcmc(G.N, T, n_spikes, ref_duration);

                % Force microstate to uint8 for snapshots (and in general)
                if ~isa(microstate0, 'uint8')
                    microstate0 = uint8(microstate0);
                end

                if do_snap
                    [x, Efin, snaps] = mcmc_core_xtsnap(microstate0, G, B, beta, n_iter, ref_duration, phi, n_snap);
                else
                    [x, Efin] = mcmc_core_xtsnap(microstate0, G, B, beta, n_iter, ref_duration, phi);
                    snaps = [];
                end

                % quick stats for the first few cells
                if (ib <= 1) && (is <= 2)
                    fprintf(1, "xstats(min=%g max=%g mean=%g) ", min(x), max(x), mean(x));
                end

                % Sanity checks
                if numel(x) ~= T
                    bad_len_x(ib, is) = true;
                    exit_code(ib, is) = int32(1001);
                    fprintf(1, 'FAIL bad_len_x numel=%d\n', numel(x));
                    continue
                end
                if any(~isfinite(x))
                    nan_in_x(ib, is) = true;
                    exit_code(ib, is) = int32(1002);
                    fprintf(1, 'FAIL nonfinite in x\n');
                    continue
                end
                if ~isfinite(Efin)
                    exit_code(ib, is) = int32(1003);
                    fprintf(1, 'FAIL nonfinite Efin\n');
                    continue
                end

                % Store (RAM)
                x_grid(ib, is, 1:T) = reshape(single(x), 1, 1, T);
                E_final(ib, is) = Efin;

                % Store snapshots (RAM)
                if do_snap
                    % snaps expected: [ks N T] same type as microstate (uint8)
                    if ~isempty(snaps)
                        if ~isa(snaps, 'uint8')
                            snaps = uint8(snaps);
                        end
                        ks = size(snaps, 1);
                        if ks > max_snaps
                            ks = max_snaps; % safety
                        end
                        micro_snaps_grid(ib, is, 1:ks, :, :) = snaps(1:ks, :, :);
                        snap_count(ib, is) = int32(ks);
                    else
                        snap_count(ib, is) = int32(0);
                    end
                end

                ok(ib, is) = true;
                exit_code(ib, is) = int32(0);

                if do_snap
                    fprintf(1, 'OK E=%g snaps=%d\n', Efin, snap_count(ib,is));
                else
                    fprintf(1, 'OK E=%g\n', Efin);
                end

            catch ME
                ok(ib, is) = false;
                exit_code(ib, is) = int32(2000);

                fprintf(2, 'EXCEPTION ib=%d is=%d beta=%g sdel=%g | %s | %s\n', ...
                        ib, is, beta, delsd, ME.identifier, ME.message);

                for kk = 1:numel(ME.stack)
                    st = ME.stack(kk);
                    fprintf(2, '  at %s (line %d)\n', st.name, st.line);
                end

                fprintf(1, 'EXCEPTION\n');
            end
        end

        fprintf(1, "--- finished beta index %d/%d ---\n", ib, nb);
    end

    % ---- hard diagnostics right before save ----
    ok_sum = nnz(ok);
    E_finite = nnz(isfinite(E_final));
    x_nnz = nnz(x_grid);
    fprintf(1, "ABOUT_TO_SAVE: ok_sum=%d/%d E_finite=%d/%d x_nnz=%d/%d\n", ...
            ok_sum, numel(ok), E_finite, numel(E_final), x_nnz, numel(x_grid));

    if do_snap
        snap_total = sum(snap_count(:));
        fprintf(1, "SNAP_DIAG: total_snaps=%d max_snaps_per_cell=%d\n", snap_total, max_snaps);
    end

    fprintf(1, "E_final finite min/max (if any): ");
    if E_finite > 0
        fprintf(1, "%g / %g\n", min(E_final(isfinite(E_final))), max(E_final(isfinite(E_final))));
    else
        fprintf(1, "(none)\n");
    end

    % Optional: text backup for E_final
    try
        txtname = sprintf('%s_Efinal.txt', prefix);
        fid = fopen(txtname, 'w');
        if fid ~= -1
            fprintf(fid, "# ib is E_final\n");
            for ib = 1:nb
                for is = 1:ns
                    fprintf(fid, "%d %d %.17g\n", ib, is, E_final(ib,is));
                end
            end
            fclose(fid);
            fprintf(1, "WROTE %s\n", txtname);
        else
            fprintf(1, "WARN: cannot open %s for writing\n", txtname);
        end
    catch ME
        fprintf(2, "WARN: Efinal txt write failed: %s | %s\n", ME.identifier, ME.message);
    end

    % ---- single-shot save (MCR-safe) ----
    fprintf(1, "ABOUT_TO_SAVE_MAT %s\n", out_file);
    save(out_file, ...
        'betas','sdels','T','N','p_conn','delmean','n_spikes','ref_duration', ...
        'n_iter','n_snap','phi_name','seed', ...
        'x_grid','E_final','ok','exit_code','nan_in_x','bad_len_x', ...
        'micro_snaps_grid','snap_count', ...
        'G_s_ind_grid','G_delay_vec_grid','G_E_grid','G_D_grid', ...
        'code_version','code_pwd','code_when','slurm_job_id', ...
        '-v7.3');

    fprintf(1, "SAVE_DONE %s\n", out_file);

end

function x = str2double_if_needed(v)
    if ischar(v) || isstring(v)
        x = str2double(v);
    else
        x = v;
    end
end
