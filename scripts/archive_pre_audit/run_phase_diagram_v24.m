function run_phase_diagram(varargin)
% run_phase_diagram PREFIX N P_CONN DELMEAN SDEL_MIN SDEL_MAX SDEL_STEP T N_SPIKES REF_DURATION ...
%                  BETA_MIN BETA_MAX BETA_STEP N_ITER PHI [SEED]
%
% Example:
% run_phase_diagram('G1',100,0.2,5, 0.0,2.0,0.2, 200,200,3, 0.2,2.0,0.2, 20000,'tanh',123)

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
    if nargin >= 15, phi_name = char(varargin{15}); end
    if nargin >= 16, seed = str2double_if_needed(varargin{16}); end

    if ~isempty(seed) && isfinite(seed)
        rng(seed, 'twister');
    end

    % ---- grids ----
    betas = 10.^(beta_min:beta_step:beta_max);
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
    code_version = 'run_phase_diagram DEBUG 2026-02-04 RAM_SAVE_V73';
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

    % IMPORTANT: for debugging save REAL values, not uint16 that might floor to zero
    x_grid = zeros(nb, ns, T, 'single');   % <--- was uint16; keep single for now
    E_final = NaN(nb, ns);

    % ---- main loop ----
    fprintf(1, "SAVE_MODE = RAM_THEN_SAVE_V73\n");
    fprintf(1, "START %s pwd=%s out_file=%s nb=%d ns=%d T=%d SLURM_JOB_ID=%s\n", ...
            datestr(now), pwd, out_file, nb, ns, T, slurm_job_id);

    for ib = 1:nb
        beta = betas(ib);
        fprintf(1, "--- ENTER ib=%d/%d beta=%g ---\n", ib, nb, beta);

        for is = 1:ns
            delsd = sdels(is);

            fprintf(1, 'ib=%d/%d beta=%g  is=%d/%d sdel=%g ... ', ...
                    ib, nb, beta, is, ns, delsd);

            try
                G = initialize_network(N, p_conn, delmean, delsd);
                B = transfer_operator(G.Adj, G.delays, T);
                microstate0 = initialize_mcmc(G.N, T, n_spikes, ref_duration);

                [x, Efin] = mcmc_core_xt(microstate0, G, B, beta, n_iter, ref_duration, phi);

                % quick stats for the first few cells (to catch silent "all <1 then uint16->0")
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
                    % treat non-finite energy as failure too
                    exit_code(ib, is) = int32(1003);
                    fprintf(1, 'FAIL nonfinite Efin\n');
                    continue
                end

                % Store (RAM)
                x_grid(ib, is, 1:T) = reshape(single(x), 1, 1, T);
                E_final(ib, is) = Efin;
                ok(ib, is) = true;
                exit_code(ib, is) = int32(0);

                fprintf(1, 'OK E=%g\n', Efin);

            catch ME
                ok(ib, is) = false;
                exit_code(ib, is) = int32(2000);

                fprintf(2, 'EXCEPTION ib=%d is=%d beta=%g sdel=%g | %s | %s\n', ...
                        ib, is, beta, delsd, ME.identifier, ME.message);

                for k = 1:numel(ME.stack)
                    st = ME.stack(k);
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
    x_nnz = nnz(x_grid);  % counts nonzero elements in entire 3D array
    fprintf(1, "ABOUT_TO_SAVE: ok_sum=%d/%d E_finite=%d/%d x_nnz=%d/%d\n", ...
            ok_sum, numel(ok), E_finite, numel(E_final), x_nnz, numel(x_grid));
    fprintf(1, "E_final finite min/max (if any): ");
    if E_finite > 0
        fprintf(1, "%g / %g\n", min(E_final(isfinite(E_final))), max(E_final(isfinite(E_final))));
    else
        fprintf(1, "(none)\n");
    end

    % Optional: text backup for E_final (useful if MAT save acts up)
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
        'n_iter','phi_name','seed', ...
        'x_grid','E_final','ok','exit_code','nan_in_x','bad_len_x', ...
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
