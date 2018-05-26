% [INPUT]
% data = A numeric array representing the sample on which the Benford's Law analysis must be performed.
% a    = A float [0.01,0.10] representing the statistical significance threshold for the tests (optional, default=0.05).

function benford_analyse(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    a = res.a;
    
    benford_analyse_internal(data,a);

end

function benford_analyse_internal(data,a)


    bd2 = benford_data(data,1);
    [mant_test,mant_desc] = benford_mantissae(bd2,a);
    [res_fo,res_so] = all_gofs(bd2,a);
    
    z = norminv(1 - (a / 2));
    s = bd2.FirstOrderTable.TheP .* (1 - bd2.FirstOrderTable.TheP);
    k = s ./ 9;
    
    intval = mean(bd2.FirstOrderTable.TheP) - (z .* (std(bd2.FirstOrderTable.TheP) ./ sqrt(1:9).'));

end

function [res_fo,res_so] = all_gofs(bd,a)

    gofs = {'AD' 'CV' 'DC' 'DE' 'FR' 'G2' 'J2' 'JD' 'JS' 'KS' 'KU' 'T2' 'U2' 'X2'};
    gofs_len = numel(gofs);
    
    res_fo = cell2table(cell(gofs_len,3),'RowNames',gofs,'VariableNames',{'H0' 'Stat' 'PValue'});
    res_so = cell2table(cell(gofs_len,3),'RowNames',gofs,'VariableNames',{'H0' 'Stat' 'PValue'});
    res_off = 1;
    
    for i = 1:gofs_len
        [h0,stat,pval] = benford_gof(bd,gofs{i},a,false);
        res_fo{res_off,:} = {h0 stat pval};
        
        [h0,stat,pval] = benford_gof(bd,gofs{i},a,true);
        res_so{res_off,:} = {h0 stat pval};

        res_off = res_off + 1;
    end
    
    res_fo.H0 = logical(res_fo.H0);
    res_so.H0 = logical(res_so.H0);

end