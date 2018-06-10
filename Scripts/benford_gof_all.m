% [INPUT]
% tab  = A table returned by the "benford_digits" function.
% a    = A float [0.01,0.10] representing the statistical significance threshold for the test (optional, default=0.05).
% sims = An integer representing the number of Monte Carlo simulations to perform (optional, default=10000).
%
% [OUTPUT]
% gofs = A 14-by-3 table containing the goodness-of-fit test results, with the following columns:
%         > H0: booleans indicating whether the null hypothesis of each test is accepted (true) or rejected (false).
%         > Statistic: the value of each statistic.
%         > pValue: the p-value associated to each statistic.

function gofs = benford_gof_all(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('tab',@(x)validateattributes(x,{'table'},{}));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
        p.addOptional('sims',10000,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',1000}));
    end

    p.parse(varargin{:});

    res = p.Results;
    tab = res.tab;
    a = res.a;
    sims = res.sims;

    gofs = benford_gof_all_internal(tab,a,sims);

end

function gofs = benford_gof_all_internal(tab,a,sims)

    gofs_tab = {
        'AD' false;
        'CV' false;
        'DC' true;
        'DE' true;
        'FR' true;
        'G2' false;
        'J2' true;
        'JD' false;
        'JS' false;
        'KS' false;
        'KU' false;
        'T2' false;
        'U2' false;
        'X2' false
    };

    gofs_len = size(gofs_tab,1);
    h0 = false(gofs_len,1);
    stat = NaN(gofs_len,1);
    pval = NaN(gofs_len,1);

    off = 1;
    
    for i = 1:gofs_len
        if (gofs_tab{i,2})
            [h0_curr,stat_curr,pval_curr] = benford_gof(tab,gofs_tab{i,1},a,sims);
        else
            [h0_curr,stat_curr,pval_curr] = benford_gof(tab,gofs_tab{i,1},a);
        end
        
        h0(off) = h0_curr;
        stat(off) = stat_curr;
        pval(off) = pval_curr;

        off = off + 1;
    end
    
    gofs = table(h0,stat,pval);
    gofs.Properties.RowNames = gofs_tab(:,1);
    gofs.Properties.VariableNames = {'H0' 'Statistics' 'pValues'};

end