repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(repo_root);

out_dir = fullfile(repo_root, 'figures', 'phase_entropy_h3');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

files = {
    fullfile(repo_root, 'results', 'raw', 'beta_0_2', ...
        'Gphase30x_betalin0to2_h3_phase.mat')
    fullfile(repo_root, 'results', 'raw', 'beta_2p2_4', ...
        'Gphase30x_betalin2p2to4_h3_phase.mat')
    fullfile(repo_root, 'results', 'raw', 'beta_4p2_6', ...
        'Gphase30x_betalin4p2to6_h3_phase.mat')
};

all_betas = [];
all_sdels = [];
all_H = [];
all_IPR = [];

for k = 1:numel(files)
    fprintf('Loading %s\n', files{k});
    S = load(files{k});

    betas = S.betas(:);
    sdels = S.sdels(:).';
    x_grid = double(S.x_grid);
    ok = logical(S.ok);

    nb = numel(betas);
    ns = numel(sdels);
    T = size(x_grid, 3);

    Hnorm = NaN(nb, ns);
    IPR = NaN(nb, ns);

    for ib = 1:nb
        for is = 1:ns
            if ~ok(ib,is)
                continue
            end

            x = squeeze(x_grid(ib,is,:));
            x = x(:);

            sx = sum(x);
            if sx <= 0
                continue
            end

            p = x ./ sx;
            ppos = p(p > 0);

            Hnorm(ib,is) = -sum(ppos .* log(ppos)) ./ log(T);
            IPR(ib,is) = sum(p .* p);
        end
    end

    all_betas = [all_betas; betas];

    if isempty(all_sdels)
        all_sdels = sdels;
    else
        if any(abs(all_sdels - sdels) > 1e-12)
            error('sdel grids do not match.');
        end
    end

    all_H = [all_H; Hnorm];
    all_IPR = [all_IPR; IPR];
end

[all_betas, idx] = sort(all_betas);
all_H = all_H(idx,:);
all_IPR = all_IPR(idx,:);

save(fullfile(out_dir, 'entropy_ipr_processed.mat'), ...
    'all_betas', 'all_sdels', 'all_H', 'all_IPR');

% ---- entropy heatmap ----
fig = figure('Color','w','Position',[100 100 900 650]);
imagesc(all_betas, all_sdels, all_H.');
set(gca,'YDir','normal');
xlabel('\beta');
ylabel('delay std s_{\Delta}');
title('Normalized temporal entropy H/log(T)');
cb = colorbar;
ylabel(cb, 'H/log(T)');
axis tight;
exportgraphics(fig, fullfile(out_dir, 'heatmap_entropy_normalized.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(out_dir, 'heatmap_entropy_normalized.pdf'), 'ContentType', 'vector');
close(fig);

% ---- IPR heatmap ----
fig = figure('Color','w','Position',[100 100 900 650]);
imagesc(all_betas, all_sdels, all_IPR.');
set(gca,'YDir','normal');
xlabel('\beta');
ylabel('delay std s_{\Delta}');
title('Temporal concentration IPR');
cb = colorbar;
ylabel(cb, 'IPR');
axis tight;
exportgraphics(fig, fullfile(out_dir, 'heatmap_ipr.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(out_dir, 'heatmap_ipr.pdf'), 'ContentType', 'vector');
close(fig);

fprintf('Saved figures to %s\n', out_dir);
