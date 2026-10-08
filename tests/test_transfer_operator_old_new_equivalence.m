function test_transfer_operator_old_new_equivalence()
%TEST_TRANSFER_OPERATOR_OLD_NEW_EQUIVALENCE
% Check that the audited explicit causal-delay operator is equivalent to
% the old implementation that relied on MATLAB's conjugating transpose '.

    repo_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(repo_root, 'src')));

    N = 5;
    K = 16;

    rng(1);

    A = rand(N) < 0.3;
    A(1:N+1:end) = 0;

    D = randi([1,5], N, N);
    D(~A) = 0;

    tt = 0:K-1;
    kern = fft(exp(-tt./5) - 0.5 .* exp(-(K-1-tt)./5));

    B_old = cell(1,K);
    B_new = cell(1,K);

    for k = 1:K
        q = k - 1;

        % Old implementation:
        % plus phase inside, then conjugating transpose ' flips the sign.
        B_old{k} = sparse(kern(k) .* ...
            (A .* exp( 2 .* pi .* 1i .* q ./ K .* D))' );

        % New implementation:
        % explicit causal negative phase and non-conjugating transpose .'
        B_new{k} = sparse(kern(k) .* ...
            (A .* exp(-2 .* pi .* 1i .* q ./ K .* D)).' );
    end

    max_err = 0;
    for k = 1:K
        err = norm(full(B_old{k} - B_new{k}), 'fro');
        max_err = max(max_err, err);
    end

    fprintf('max old-new transfer operator error: %.3e\n', max_err);

    if max_err > 1e-12
        error('test_transfer_operator_old_new_equivalence: old and new operators differ.');
    end

    fprintf('test_transfer_operator_old_new_equivalence OK\n');

end
