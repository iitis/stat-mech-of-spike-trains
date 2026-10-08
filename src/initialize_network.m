function G = initialize_network(N, p_conn, delmean, delsd, seed)
%
% G = initialize_network(N, p_conn, delmean, delsd [, seed])

    if nargin < 5
        seed = [];
    end
    if ~isempty(seed) && isfinite(seed)
        rng(seed, 'twister');
    end

    % --- graph ---
    G.name = sprintf('mcmc_net_(%g,%g,%g).mat', p_conn, delmean, delsd);
    G.N = N;

    Adj = rand(N) < p_conn;
    Adj(1:N+1:end) = 0;         % remove self-connections
    G.Adj = Adj ~= 0;
    G.E = nnz(G.Adj);

    % --- edges ---
    G.s_ind = find(G.Adj);
    m = numel(G.s_ind);

    % --- delays ---
    if m == 0
        G.delays = zeros(N, N);
        G.D = 0;
        fprintf('Network: N=%d p=%g E=0 (no edges)\n', N, p_conn);
        return;
    end

    if delsd < 0.5 * delmean
        sample = rounded_randn(delmean, delsd, m);  % as column
        if any(sample == 1)
            sample = ceil(me_sample(delmean - 0.5, delsd, m));
        else
            disp('norm*');
        end
    else
        sample = ceil(me_sample(delmean - 0.5, delsd, m));
        if any(~isfinite(sample))
            sample = ceil(exprnd(delmean - 0.5, m, 1));
            disp('exp*');
        end
    end

    sample = sample(:);
    sample = max(sample, 1);  % delay >= 1

    fprintf('mean delay = %g   delay std = %g\n', mean(sample), std(sample));

    G.delays = zeros(N, N);
    G.delays(G.s_ind) = sample;
    G.D = full(max(G.delays(:)));
end
