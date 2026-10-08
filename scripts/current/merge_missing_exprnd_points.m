repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(repo_root);

specs = {
    fullfile(repo_root,'results','raw','beta_0_2','Gphase30x_betalin0to2_h3_phase.mat'), ...
    fullfile(repo_root,'results','raw','missing_exprnd_fix','b0_beta0p2_sdel4p5_phase.mat'), 0.2, 4.5

    fullfile(repo_root,'results','raw','beta_0_2','Gphase30x_betalin0to2_h3_phase.mat'), ...
    fullfile(repo_root,'results','raw','missing_exprnd_fix','b0_beta1p8_sdel4p5_phase.mat'), 1.8, 4.5

    fullfile(repo_root,'results','raw','beta_0_2','Gphase30x_betalin0to2_h3_phase.mat'), ...
    fullfile(repo_root,'results','raw','missing_exprnd_fix','b0_beta2p0_sdel4p5_phase.mat'), 2.0, 4.5

    fullfile(repo_root,'results','raw','beta_2p2_4','Gphase30x_betalin2p2to4_h3_phase.mat'), ...
    fullfile(repo_root,'results','raw','missing_exprnd_fix','b2_beta2p4_sdel4p5_phase.mat'), 2.4, 4.5

    fullfile(repo_root,'results','raw','beta_2p2_4','Gphase30x_betalin2p2to4_h3_phase.mat'), ...
    fullfile(repo_root,'results','raw','missing_exprnd_fix','b2_beta4p0_sdel4p5_phase.mat'), 4.0, 4.5

    fullfile(repo_root,'results','raw','beta_4p2_6','Gphase30x_betalin4p2to6_h3_phase.mat'), ...
    fullfile(repo_root,'results','raw','missing_exprnd_fix','b4_beta4p4_sdel4p5_phase.mat'), 4.4, 4.5

    fullfile(repo_root,'results','raw','beta_4p2_6','Gphase30x_betalin4p2to6_h3_phase.mat'), ...
    fullfile(repo_root,'results','raw','missing_exprnd_fix','b4_beta6p0_sdel4p5_phase.mat'), 6.0, 4.5
};

for k = 1:size(specs,1)
    base_file  = specs{k,1};
    patch_file = specs{k,2};
    beta_val   = specs{k,3};
    sdel_val   = specs{k,4};

    fprintf('\n=== MERGE %d/%d ===\n', k, size(specs,1));
    fprintf('base:  %s\n', base_file);
    fprintf('patch: %s\n', patch_file);
    fprintf('beta=%g sdel=%g\n', beta_val, sdel_val);

    B = load(base_file);
    P = load(patch_file);

    ib = find(abs(B.betas - beta_val) < 1e-10, 1);
    is = find(abs(B.sdels - sdel_val) < 1e-10, 1);

    if isempty(ib) || isempty(is)
        error('Could not locate beta=%g sdel=%g in base file.', beta_val, sdel_val);
    end

    if ~P.ok(1,1)
        error('Patch file is not OK: %s', patch_file);
    end

    fprintf('Target cell: ib=%d is=%d\n', ib, is);

    % ---- core numerical arrays ----
    B.x_grid(ib,is,:) = P.x_grid(1,1,:);
    B.E_final(ib,is) = P.E_final(1,1);
    B.ok(ib,is) = P.ok(1,1);
    B.exit_code(ib,is) = P.exit_code(1,1);
    B.nan_in_x(ib,is) = P.nan_in_x(1,1);
    B.bad_len_x(ib,is) = P.bad_len_x(1,1);

    % ---- diagnostics ----
    B.final_spike_count_grid(ib,is) = P.final_spike_count_grid(1,1);
    B.spike_count_ok_grid(ib,is) = P.spike_count_ok_grid(1,1);
    B.refractory_ok_grid(ib,is) = P.refractory_ok_grid(1,1);
    B.accept_rate_grid(ib,is) = P.accept_rate_grid(1,1);
    B.n_accept_grid(ib,is) = P.n_accept_grid(1,1);
    B.n_reject_grid(ib,is) = P.n_reject_grid(1,1);

    % ---- snapshots ----
    B.micro_snaps_grid(ib,is,:,:,:) = P.micro_snaps_grid(1,1,:,:,:);
    B.snap_count(ib,is) = P.snap_count(1,1);

    % ---- optional graph fields ----
    graph_fields = {'G_s_ind_grid','G_delay_vec_grid','G_E_grid','G_D_grid'};

    for gf = 1:numel(graph_fields)
        fname = graph_fields{gf};

        if isfield(B, fname) && isfield(P, fname)
            B.(fname) = assign_patch_value(B.(fname), P.(fname), ib, is, fname);
        else
            fprintf('Skipping missing optional field: %s\n', fname);
        end
    end

    % ---- patch history ----
    if ~isfield(B, 'patch_history') || ~iscell(B.patch_history)
        B.patch_history = {};
    end

    B.patch_history{end+1,1} = struct( ...
        'patch_file', patch_file, ...
        'beta', beta_val, ...
        'sdel', sdel_val, ...
        'merged_at', char(datetime('now')));

    fprintf('Saving patched base file...\n');
    save(base_file, '-struct', 'B', '-v7.3');
    fprintf('MERGE_DONE beta=%g sdel=%g into ib=%d is=%d\n', beta_val, sdel_val, ib, is);
end

fprintf('\nAll missing points merged.\n');

function base_field = assign_patch_value(base_field, patch_field, ib, is, fname)
%ASSIGN_PATCH_VALUE Robustly assign one-point patch into a grid field.
%
% Handles both cell arrays and numeric / matrix arrays. For one-point patch
% files, some fields may not be cell arrays. Because MATLAB enjoys making
% adults sad, we handle both cases explicitly.

    if iscell(base_field)
        if iscell(patch_field)
            val = patch_field{1,1};
        else
            val = patch_field;
        end

        base_field{ib,is} = val;
        fprintf('Merged cell field %s\n', fname);
        return
    end

    % Non-cell base field.
    if iscell(patch_field)
        val = patch_field{1,1};
    else
        val = patch_field;
    end

    szb = size(base_field);
    szv = size(val);

    try
        if isscalar(val)
            base_field(ib,is) = val;
        elseif ndims(base_field) == 2 && isequal(size(base_field), [numel(szb) numel(szb)])
            base_field(ib,is) = val;
        elseif ndims(base_field) >= 3
            % Generic attempt for arrays indexed by (ib,is,...).
            idx = repmat({':'}, 1, ndims(base_field));
            idx{1} = ib;
            idx{2} = is;
            base_field(idx{:}) = val;
        else
            % If this is a graph-like matrix field rather than an indexed grid,
            % do not force a dangerous assignment.
            fprintf('WARN: skipped non-cell field %s with base size [%s] and patch size [%s]\n', ...
                fname, num2str(szb), num2str(szv));
        end
        fprintf('Merged non-cell field %s\n', fname);
    catch ME
        fprintf('WARN: could not merge optional field %s: %s | %s\n', ...
            fname, ME.identifier, ME.message);
    end
end
