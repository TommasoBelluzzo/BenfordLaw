% [INPUT]
% data = A numeric array of n elements representing the sample on which the Duplication Analysis must be performed.
%
% [OUTPUT]
% df   = A 1-by-5 table containing the Distortion Factor test result, with the following columns:
%         > AM: a float representing the actual mean of the collapsed values.
%         > EM: a float representing the expected mean of the collapsed values.
%         > Value: a float representing the distortion factor.
%         > pValue: a float representing the p-value associated to the distortion factor.
%         > Significance: a boolean indicating whether the test is significant (true) or not (false) for the given confidence level.
% coll = A numeric vector representing the collapsed values of the sample.

function [d10n,d0n,d0p,d10p] = benford_duplication(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    
    if (nargout ~= 4)
        error('Only 4 output arguments can be specified.');
    end

    [d10n,d0n,d0p,d10p] = benford_duplication_internal(data);

end

function [d10n,d0n,d0p,d10p] = benford_duplication_internal(data)

    data_10n = data(data <= -10);

    if (numel(data_10n) > 0)
        data_10n_uni = unique(data_10n);
        [a,b] = histcounts(data_10n,data_10n_uni);
    end

    data_0n = data((data > -10) & (data < 0));
    
    if (numel(data_0n) > 0)
        data_0n_uni = unique(data_0n);
        [a,b] = histcounts(data_0n,data_0n_uni);
    end

    data_0p = data((data > 0) & (data < 10));
    
    if (numel(data_0p) > 0)
        data_0p = sort(data_0p);

    end

    data_10p = data(data >= 10);
    
    if (numel(data_10p) > 0)
        [unique_values, ~, labels] = unique(data_10p);
        counts = accumarray(labels, 1);
    end

end