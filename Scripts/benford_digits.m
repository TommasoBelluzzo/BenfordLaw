% [INPUT]
% data = A numeric array of n elements representing the sample on which the specified digits extraction must be performed.
% extr = A string representing the digits extraction to perform, its value can be one of the following:
%         - 1ST (first digit extraction)
%         - 2ND (second digit extraction)
%         - 3RD (third digit extraction)
%         - F2D (first-two digits extraction)
%         - F3D (first-three digits extraction)
%         - L2D (last-two digits extraction)
% a    = A float [0.01,0.10] representing the statistical significance threshold for the Z-score tests (optional, default=0.05).
% ccf  = A boolean indicating whether to apply a continuity correction factor on Z-scores (optional, default=true).
% ran  = A string representing the range of values to consider (optional, default='ALL').
%        Its value can be one of the following:
%         - ALL (all values)
%         - NEG (only negative values)
%         - POS (only positive values)
% dec  = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
% btv  = A boolean indicating whether to perform data transformation and validation (optional, default=true).
%        This input argument should be set to false only when data has already been transformed and validated by another function.
%
% [OUTPUT]
% tab = 

function tab = benford_digits(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addRequired('extr',@(x)any(validatestring(x,{'1ST','2ND','3RD','F2D','F3D','L2D'})));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
        p.addOptional('ccf',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
        p.addOptional('ran','ALL',@(x)any(validatestring(x,{'ALL','NEG','POS'})));
        p.addOptional('dec',2,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',10}));
        p.addOptional('btv',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    extr = res.extr;
    a = res.a;
    ccf = res.ccf;
    ran = res.ran;
    dec = res.dec;
    btv = res.btv;

    if (btv)
        data = benford_data(data,ran,dec);
    end

    tab = benford_digits_internal(data,extr,a,ccf);

end

function tab = benford_digits_internal(data,extr,a,ccf)

    switch (extr)
        case 'L2D'
            [the_p,the_f,emp_c,emp_p,emp_f] = extract_last_two(data);
    end

    n = numel(the_p);

    diff_p = emp_p - the_p;
    diff_p_abs = abs(diff_p);

    z_targ = norminv(1 - (a / 2));
    z_num = diff_p_abs;
    z_den = sqrt((the_p .* (1 - the_p)) ./ n);

    the_p_lb = the_p - (z_targ .* z_den);
    the_p_ub = the_p + (z_targ .* z_den);
    
    if (ccf)
        ccf_val = 1 / (2 * n);

        ccf_idx = z_num > ccf_val;
        z_num(ccf_idx) = z_num(ccf_idx) - ccf_val;
        
        ccf_idx = the_p_lb > ccf_val;
        the_p_lb(ccf_idx) = the_p_lb(ccf_idx) - ccf_val;

        the_p_ub = the_p_ub - ccf_val;
    end

    z = z_num ./ z_den;
    z_test = z > z_targ;

    
    tab = table(the_p_lb,the_p,the_p_ub,the_f,emp_c,emp_p,emp_f,z,z_test);
end

function [the_p,the_f,emp_c,emp_p,emp_f] = extract_last_two(data)

    the_p = repmat(1/100,100,1);
    the_f = [cumsum(the_p(1:end-1)); 1];

    d = mod(data,100);
    n = numel(d);

    emp_c = histcounts(d,[(0:99).'; Inf]).';
    emp_p = emp_c ./ n;
    emp_f = [cumsum(emp_p(1:end-1)); 1];

end