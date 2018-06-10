% [INPUT]
% data = A numeric array of n elements representing the sample on which the Distortion Factor Model must be applied.
% ddec = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
%        No rounding is performed, the exceeding decimals are truncated as if they were not present.
% a    = A float [0.01,0.10] representing the statistical significance threshold for the Mantissae Arc test (optional, default=0.05).
%
% [OUTPUT]
% df   = A 1-by-5 table containing the Distortion Factor test result, with the following columns:
%         > AM: a float representing the actual mean of the collapsed values.
%         > EM: a float representing the expected mean of the collapsed values.
%         > Value: a float representing the distortion factor.
%         > pValue: a float representing the p-value associated to the distortion factor.
%         > Significance: a boolean indicating whether the test is significant (true) or not (false) for the given confidence level.
% coll = A numeric vector representing the collapsed values of the sample.

function [df,coll] = benford_df(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('ddec',2,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',10}));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    ddec = res.ddec;
    a = res.a;

    switch (nargout)
        case 1
            df = benford_df_internal(data,ddec,a);

        case 2
            [df,coll] = benford_df_internal(data,ddec,a);
            
        otherwise
            error('Only up to 2 output arguments can be specified.');
    end

end

function [df,coll] = benford_df_internal(data,ddec,a)

    data = double(data(:));
    data = floor(data .* (10 ^ ddec)) ./ (10 ^ ddec);
    data = data(data >= 10);

    n = numel(data);

    lt = log10(data);
    lt = lt - rem(lt,1);
    
    coll_int = (10 .* data) ./ (10 .^ lt);
    
    am = mean(coll_int);
    em = 90 / (n * ((10 ^ (1 / n)) - 1));
    df_val = (am - em) / em;

    pval = 2 * (1 - normcdf(df_val / (0.638253 / sqrt(n))));
    pval = max([0 min([pval 1])]);
    
    sig = pval >= a;

    df = table(am,em,df_val,pval,sig);
    df.Properties.VariableNames = {'AM' 'EM' 'Value' 'pValue' 'Significance'};

    if (nargout == 2)
        coll = coll_int;
    end

end