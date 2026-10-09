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

out_dir = fullfile(repo_root, 'results', 'processed', 'main_map');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
base_files = unique(specs(:,1), 'stable');
fields = {'x_grid','E_final','ok','exit_code','nan_in_x','bad_len_x', ...
    'final_spike_count_grid','spike_count_ok_grid','refractory_ok_grid', ...
    'accept_rate_grid','n_accept_grid','n_reject_grid', ...
    'micro_snaps_grid','snap_count', ...
    'G_s_ind_grid','G_delay_vec_grid','G_E_grid','G_D_grid'};
parameters = {'N','T','p_conn','delmean','n_spikes','ref_duration', ...
    'n_iter','n_snap','phi_name'};

% Always reconstruct from raw inputs, never from a previous processed output.
for batch = 1:numel(base_files)
    base_file = base_files{batch};
    B = load(base_file);
    B.source_file_grid = repmat({relative_path(base_file, repo_root)}, size(B.ok));
    B.seed_grid = repmat({B.seed}, size(B.ok));
    B.reconstruction_history = {};
    rows = find(strcmp(specs(:,1), base_file));
    for k = rows(:).'
        patch_file = specs{k,2};
        beta_val = specs{k,3};
        sdel_val = specs{k,4};
        P = load(patch_file);
        ib = find(abs(B.betas - beta_val) < 1e-10);
        is = find(abs(B.sdels - sdel_val) < 1e-10);
        assert(numel(ib) == 1 && numel(is) == 1, 'Ambiguous or missing target point.');
        assert(numel(P.betas) == 1 && numel(P.sdels) == 1 && ...
            abs(P.betas - beta_val) < 1e-10 && abs(P.sdels - sdel_val) < 1e-10, ...
            'Patch coordinates do not match the target.');
        for j = 1:numel(parameters)
            name = parameters{j};
            assert(isequaln(B.(name), P.(name)), 'Patch parameter mismatch: %s', name);
        end
        assert(isscalar(P.ok) && P.ok && P.spike_count_ok_grid && ...
            P.refractory_ok_grid && P.final_spike_count_grid == B.n_spikes && ...
            P.exit_code == 0 && ~P.nan_in_x && ~P.bad_len_x && ...
            isfinite(P.E_final) && isfinite(P.accept_rate_grid) && ...
            numel(P.x_grid) == B.T && all(isfinite(P.x_grid(:))), ...
            'Invalid replacement point: %s', patch_file);
        for j = 1:numel(fields)
            name = fields{j};
            assert(isfield(B, name) && isfield(P, name), 'Missing field: %s', name);
            B.(name) = assign_patch_value(B.(name), P.(name), ib, is, name);
        end
        source = relative_path(patch_file, repo_root);
        B.source_file_grid{ib,is} = source;
        B.seed_grid{ib,is} = P.seed;
        B.reconstruction_history{end+1,1} = struct( ...
            'base_file', relative_path(base_file, repo_root), ...
            'patch_file', source, 'beta', beta_val, 'sdel', sdel_val, ...
            'seed', P.seed, 'code_version', P.code_version);
    end
    assert(all(B.ok(:)) && all(B.spike_count_ok_grid(:)) && ...
        all(B.refractory_ok_grid(:)) && all(isfinite(B.E_final(:))) && ...
        all(isfinite(B.x_grid(:))) && all(isfinite(B.accept_rate_grid(:))) && ...
        all(B.final_spike_count_grid(:) == B.n_spikes) && ...
        all(B.exit_code(:) == 0) && ~any(B.nan_in_x(:)) && ~any(B.bad_len_x(:)), ...
        'Reconstructed batch still contains invalid points.');
    [~, filename] = fileparts(base_file);
    out_file = fullfile(out_dir, [filename '.mat']);
    B.reconstruction_version = 'main_map_reconstruction_v1';
    save(out_file, '-struct', 'B', '-v7.3');
    txt_file = fullfile(out_dir, [strrep(filename, '_phase', '') '_Efinal.txt']);
    write_energy_export(txt_file, B);
    fprintf('Saved reconstructed MAT and TXT: %s\n', out_file);
end
fprintf('Main-map reconstruction complete; raw inputs were preserved.\n');

function result = relative_path(filename, root)
    result = strrep(filename(numel(root)+2:end), filesep, '/');
end

function write_energy_export(filename, B)
    fid = fopen(filename, 'w');
    assert(fid ~= -1, 'Cannot open energy export: %s', filename);
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# ib is beta sdel E_final accept_rate ok exit_code\n');
    for ib = 1:numel(B.betas)
        for is = 1:numel(B.sdels)
            fprintf(fid, '%d %d %.17g %.17g %.17g %.17g %d %d\n', ...
                ib, is, B.betas(ib), B.sdels(is), B.E_final(ib,is), ...
                B.accept_rate_grid(ib,is), B.ok(ib,is), B.exit_code(ib,is));
        end
    end
end

function base_field = assign_patch_value(base_field, patch_field, ib, is, fname)
    if iscell(base_field)
        if iscell(patch_field)
            assert(numel(patch_field) == 1, 'Expected a one-point cell field: %s', fname);
            base_field{ib,is} = patch_field{1};
        else
            % Historical one-point files may store graph vectors directly.
            base_field{ib,is} = patch_field;
        end
    else
        if iscell(patch_field)
            assert(numel(patch_field) == 1, 'Expected a one-point cell field: %s', fname);
            patch_field = patch_field{1};
        end
        assert(size(patch_field,1) == 1 && size(patch_field,2) == 1, ...
            'Expected a one-point numeric field: %s', fname);
        dims = max(ndims(base_field), ndims(patch_field));
        for d = 3:dims
            assert(size(base_field,d) == size(patch_field,d), ...
                'Incompatible trailing dimensions: %s', fname);
        end
        idx = repmat({':'}, 1, dims);
        idx{1} = ib; idx{2} = is;
        base_field(idx{:}) = patch_field;
    end
end
