% [INPUT]
% data = A numeric array of n elements representing the sample on which the Mantissae Analysis must be performed.
% ran  = A string representing the range of values to consider (optional, default='ALL').
%        Its value can be one of the following:
%         - ALL (all values)
%         - NEG (only negative values)
%         - POS (only positive values)
% a    = A float [0.01,0.10] representing the statistical significance threshold for the Mantissae Arc test (optional, default=0.05).
%
% [OUTPUT]
% mant = An n-by-3 table containing the analyzed components, with the following columns:
%         > Mantissae: a vector of floats containing the mantissae of the sample.
%         > X: a vector of floats representing the X coordinates of the mantissae.
%         > Y: a vector of floats representing the Y coordinates of the mantissae.
% test = A 1-by-3 table containing the Mantissae Arc test result, with the following columns:
%         > H0: a boolean indicating whether the null hypothesis is accepted (true) or rejected (false).
%         > Statistic: a float representing the statistic value.
%         > pValue: a float representing the p-value associated to the statistic.
% desc = A 4-by-2 table containing the theoretical and empirical descriptive statistics, with the following columns:
%         > The: the theoretical mean, variance, skewness and kurtosis.
%         > Emp: the empirical mean, variance, skewness and kurtosis.

function [mant,test,desc] = benford_mantissae(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('ran','ALL',@(x)any(validatestring(x,{'ALL','NEG','POS'})));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    ran = res.ran;
    a = res.a;
    
    data = double(data(:));
    
    switch (ran)
        case 'NEG'
            data = abs(data(isfinite(data) & (data < 0)));
    
        case 'POS'
            data = data(isfinite(data) & (data > 0));
            
        otherwise
            data = abs(data(isfinite(data) & (data ~= 0)));
    end

    if (numel(unique(data)) < 10)
        error('The number of unique valid observations in the sample must be greater than or equal to 10.');
    end
    
    if (numel(data) < 250)
        warning('A minimum sample size of 250 valid observations is recommended in order to produce a coherent analysis.');
    end

    switch (nargout)
        case 1
            mant = benford_mantissae_internal(data,a);

        case 2
            [mant,test] = benford_mantissae_internal(data,a);
            
        case 3
            [mant,test,desc] = benford_mantissae_internal(data,a);
            
        otherwise
            error('Only up to 3 output arguments can be specified.');
    end

end

function [mant,test,desc] = benford_mantissae_internal(data,a)

    data = abs(data);
    
    l = log10(data);
    
    neg_idx = l < 0;
    l(neg_idx) = l(neg_idx) + abs(ceil(l(neg_idx))) + 1;

    l_trun = str2double(cellstr(num2str(l)));
    l_trun = l_trun - rem(l_trun,1);
    
    m = l - l_trun;
    x = cos(2 .* pi() .* m);
    y = sin(2 .* pi() .* m);
    
    mant = table(m,x,y);
    mant.Properties.VariableNames = {'Mantissae' 'X' 'Y'};

    if (nargout >= 2)
        stat = (mean(x) ^ 2) + (mean(y) ^ 2);
        pval = 1 - exp(-stat * numel(m));
        h0 = pval >= a;

        test = table(h0,stat,pval);
        test.Properties.VariableNames = {'H0' 'Statistic' 'pValue'};
    end
    
    if (nargout == 3)
        the = [0.5; 0.08333; 0; -1.2];
        emp = [mean(m); var(m); skewness(m,true); kurtosis(m,true)];
        
        desc = table(the,emp);
        desc.Properties.RowNames = {'Mean' 'Variance' 'Skewness' 'Kurtosis'};
        desc.Properties.VariableNames = {'The' 'Emp'};
    end

end