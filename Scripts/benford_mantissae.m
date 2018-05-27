% [INPUT]
% bd   = An instance of the BenfordData class produced by the "benford_data" function.
% a    = A float [0.01,0.10] representing the statistical significance threshold for the test (optional, default=0.05).
%
% [OUTPUT]
% test = A 1-by-3 table containing the Mantissae Arc test result. It has the following columns:
%         > H0: a boolean indicating whether the null hypothesis is accepted (true) or rejected (false).
%         > Statistic: a float representing the statistic value.
%         > pValue: s float representing the p-value associated to the statistic.
% desc = A 4-by-2 table containing the theoretical and empirical descriptive statistics. It has the following columns:
%         > The: the theoretical mean, variance, skewness and kurtosis.
%         > The: the empirical mean, variance, skewness and kurtosis.

function [test,desc] = benford_mantissae(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('bd',@(x)validateattributes(x,{'BenfordData'},{'scalar'}));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
    end

    p.parse(varargin{:});

    res = p.Results;
    bd = res.bd;
    a = res.a;

    switch (nargout)
        case 1
            test = benford_mantissae_internal(bd,a);

        case 2
            [test,desc] = benford_mantissae_internal(bd,a);
            
        otherwise
            error('Only 1 or 2 output arguments can be specified.');
    end

end

function [test,desc] = benford_mantissae_internal(bd,a)

    mant = bd.Mantissae;

    x = cos(2 .* pi() .* mant);
    y = sin(2 .* pi() .* mant);

    stat = (mean(x) ^ 2) + (mean(y) ^ 2);
    pval = 1 - exp(-stat * numel(mant));
    h0 = pval >= a;
    
    test = table(h0,stat,pval);
    test.Properties.VariableNames = {'H0' 'Statistic' 'pValue'};
    
    if (nargout == 2)
        desc = table([0.5; 0.08333; 0; -1.2],[mean(mant); var(mant); skewness(mant,true); kurtosis(mant,true)]);
        desc.Properties.RowNames = {'Mean' 'Variance' 'Skewness' 'Kurtosis'};
        desc.Properties.VariableNames = {'The' 'Emp'};
    end

end