% [INPUT]
% data = A numeric array of n elements representing the sample on which the Zipf's Law Analysis must be performed.
% dran = A string representing the range of values to consider (optional, default='ALL').
%        Its value can be one of the following:
%         - ALL (all values)
%         - NEG (only negative values)
%         - POS (only positive values)
% ddec = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
%        No rounding is performed, the exceeding decimals are truncated as if they were not present.
% btv  = A boolean indicating whether to perform data transformation and validation (optional, default=true).
%        This input argument should be set to false only when data has already been transformed and validated by another function.
%
% [OUTPUT]
% tab  = A k-by-6 table (where k is the number of unique elements in the sample) containing the distributions, with the following columns:
%         > Rank: a vector of integers representing the rank.
%         > Sample: a vector of numeric values representing the unique elements in the sample.
%         > RankL: a vector of floats representing the logarithm of the ranks.
%         > SampleL: a vector of floats representing logarithm of the unique elements in the sample.
%         > TheP: a vector of floats representing the theoretical proportions.
%         > EmpP: a vector of floats representing the empirical proportions.
% fit  = A 3-by-2 table containing the fitting results, with the following columns:
%         > The: intercept (A), slope (B) and the coefficient of determination (R2) of the theoretical fitting.
%         > Emp: intercept (A), slope (B) and the coefficient of determination (R2) of the empirical fitting.

function [tab,fit] = benford_zipf(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('dran','ALL',@(x)any(validatestring(x,{'ALL','NEG','POS'})));
        p.addOptional('ddec',2,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',10}));
        p.addOptional('btv',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    dran = res.dran;
    ddec = res.ddec;
    btv = res.btv;

    if (btv)
        data = benford_data(data,dran,ddec);
    end

    switch (nargout)
        case 1
            tab = benford_zipf_internal(data);

        case 2
            [tab,fit] = benford_zipf_internal(data);
            
        otherwise
            error('Only up to 2 output arguments can be specified.');
    end

end

function [tab,fit] = benford_zipf_internal(data)

	data = sort(unique(data),'descend');
    data_log = log(data);
    
	n = numel(data);

    rnk = (1:n).';
    rnk_log = log(rnk);

    emp_p = data_log ./ max(data_log);
    the_p = 1 ./ rnk;

    tab = table(rnk,data,rnk_log,data_log,the_p,emp_p);
    tab.Properties.VariableNames = {'Rank' 'Sample' 'RankL' 'SampleL' 'TheP' 'EmpP'};

    if (nargout == 2)
        the_a = data_log(1); 
        the_b = -1 * (the_a / log(n));
        the_l = the_a + (the_b .* rnk_log);

        the_rss = sum((the_l - data_log) .^ 2);
        the_tss = (n - 1) * var(the_l);
        the_r2 = min([max([(1 - (the_rss / the_tss)) 0]) 1]);

        fit = polyfit(rnk_log,data_log,1);
        emp_a = fit(2);
        emp_b = fit(1);
        emp_l = emp_a + (emp_b .* rnk_log);
        
        emp_rss = sum((data_log - emp_l) .^ 2);
        emp_tss = (n - 1) * var(emp_l);
        emp_r2 = min([max([(1 - (emp_rss / emp_tss)) 0]) 1]);
        
        fit = table([the_a; the_b; the_r2],[emp_a; emp_b; emp_r2]);
        fit.Properties.RowNames = {'A' 'B' 'R2'};
        fit.Properties.VariableNames = {'The' 'Emp'};
    end

end