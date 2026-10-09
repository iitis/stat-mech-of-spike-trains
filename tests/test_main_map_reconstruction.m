function test_main_map_reconstruction()
% Verify historical replacements, unchanged cells, provenance, and TXT export.
% Run after merge_missing_exprnd_points.m; not part of sampler-only tests.
    root = fileparts(fileparts(mfilename('fullpath')));
    batches = {'beta_0_2','beta_2p2_4','beta_4p2_6'};
    prefixes = {'Gphase30x_betalin0to2_h3', ...
        'Gphase30x_betalin2p2to4_h3','Gphase30x_betalin4p2to6_h3'};
    replacements = { ...
        {'b0_beta0p2_sdel4p5','b0_beta1p8_sdel4p5','b0_beta2p0_sdel4p5'}, ...
        {'b2_beta2p4_sdel4p5','b2_beta4p0_sdel4p5'}, ...
        {'b4_beta4p4_sdel4p5','b4_beta6p0_sdel4p5'}};
    fields = {'x_grid','E_final','ok','exit_code','nan_in_x','bad_len_x', ...
        'final_spike_count_grid','spike_count_ok_grid','refractory_ok_grid', ...
        'accept_rate_grid','n_accept_grid','n_reject_grid','snap_count', ...
        'G_s_ind_grid','G_delay_vec_grid','G_E_grid','G_D_grid'};
    for batch = 1:3
        raw_rel = ['results/raw/' batches{batch} '/' prefixes{batch} '_phase.mat'];
        out_file = fullfile(root,'results','processed','main_map', ...
            [prefixes{batch} '_phase.mat']);
        raw_file = fullfile(root, raw_rel);
        expected = load(raw_file, fields{:}, 'betas','sdels','seed');
        actual = load(out_file, fields{:}, 'betas','sdels', ...
            'source_file_grid','seed_grid','reconstruction_history');
        sources = repmat({raw_rel}, size(expected.ok));
        seeds = repmat({expected.seed}, size(expected.ok));
        out_mat = matfile(out_file);
        assert(numel(actual.reconstruction_history) == numel(replacements{batch}), ...
            'Unexpected reconstruction history length.');
        for k = 1:numel(replacements{batch})
            patch_rel = ['results/raw/missing_exprnd_fix/' replacements{batch}{k} '_phase.mat'];
            patch_file = fullfile(root, patch_rel);
            P = load(patch_file, fields{:}, 'betas','sdels','seed');
            ib = find(abs(expected.betas - P.betas) < 1e-10);
            is = find(abs(expected.sdels - P.sdels) < 1e-10);
            assert(numel(ib) == 1 && numel(is) == 1, 'Invalid test coordinates.');
            for j = 1:numel(fields)
                name = fields{j};
                value = P.(name);
                if iscell(expected.(name))
                    if iscell(value), value = value{1}; end
                    expected.(name){ib,is} = value;
                else
                    if iscell(value), value = value{1}; end
                    idx = repmat({':'}, 1, ndims(expected.(name)));
                    idx{1} = ib; idx{2} = is;
                    expected.(name)(idx{:}) = value;
                end
            end
            patch_mat = matfile(patch_file);
            assert(isequaln(out_mat.micro_snaps_grid(ib,is,:,:,:), ...
                patch_mat.micro_snaps_grid(1,1,:,:,:)), 'Replacement snapshots differ.');
            sources{ib,is} = patch_rel;
            seeds{ib,is} = P.seed;
            history = actual.reconstruction_history{k};
            assert(strcmp(history.patch_file, patch_rel) && ...
                strcmp(history.base_file, raw_rel) && isequaln(history.seed, P.seed) && ...
                abs(history.beta - P.betas) < 1e-10 && ...
                abs(history.sdel - P.sdels) < 1e-10, 'Incorrect patch history.');
        end
        for j = 1:numel(fields)
            name = fields{j};
            assert(isequaln(actual.(name), expected.(name)), ...
                'Unexpected changes or incorrect replacement in %s.', name);
        end
        assert(isequal(actual.betas,expected.betas) && ...
            isequal(actual.sdels,expected.sdels), 'Grid coordinates changed.');
        assert(isequal(actual.source_file_grid, sources) && ...
            isequaln(actual.seed_grid, seeds), 'Incorrect per-point provenance.');
        assert(all(actual.ok(:)) && all(isfinite(actual.E_final(:))), 'Invalid reconstructed map.');
        txt = readmatrix(fullfile(root,'results','processed','main_map', ...
            [prefixes{batch} '_Efinal.txt']), 'FileType','text','CommentStyle','#');
        expected_txt = zeros(numel(actual.ok),8);
        row = 0;
        for ib = 1:numel(actual.betas)
            for is = 1:numel(actual.sdels)
                row = row + 1;
                expected_txt(row,:) = [ib,is,actual.betas(ib),actual.sdels(is), ...
                    double(actual.E_final(ib,is)),double(actual.accept_rate_grid(ib,is)), ...
                    double(actual.ok(ib,is)),double(actual.exit_code(ib,is))];
            end
        end
        assert(isequaln(txt, expected_txt), 'TXT export differs from reconstructed MAT.');
        fprintf('Reconstruction verified: %s (%d points, %d replacements)\n', ...
            prefixes{batch}, numel(actual.ok), numel(replacements{batch}));
    end
    fprintf('test_main_map_reconstruction OK\n');
end
