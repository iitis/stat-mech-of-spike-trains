function r = exprnd(mu, varargin)
%EXPRND Minimal exponential random number generator.
%
% Local replacement for Statistics Toolbox exprnd.
% Generates exponential random variables with mean mu using inverse CDF:
%
%     X = -mu * log(U),   U ~ Uniform(0,1).

    if nargin < 1
        error('exprnd:NotEnoughInputs', 'Mean parameter mu is required.');
    end

    if any(mu(:) < 0)
        error('exprnd:BadMean', 'Mean parameter mu must be non-negative.');
    end

    if isempty(varargin)
        if isscalar(mu)
            sz = [1 1];
        else
            sz = size(mu);
        end
    elseif numel(varargin) == 1 && isnumeric(varargin{1})
        sz = varargin{1};
    else
        sz = cell2mat(varargin);
    end

    U = rand(sz);
    U = max(U, realmin);

    if isscalar(mu)
        r = -mu .* log(U);
    else
        if ~isequal(size(mu), sz)
            error('exprnd:SizeMismatch', ...
                'Non-scalar mu must match requested output size.');
        end
        r = -mu .* log(U);
    end
end
