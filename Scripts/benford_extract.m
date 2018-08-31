% [INPUT]
% data = A numeric array of n elements representing the sample on which the specified digits extraction must be performed.
% dran = A string representing the range of values to consider (optional, default='ALL').
%        Its value can be one of the following:
%         - ALL (all values)
%         - NEG (only negative values)
%         - POS (only positive values)
% ddec = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
%        No rounding is performed, the exceeding decimals are truncated as if they were not present.
% extr = A string representing the digits extraction to perform (optional, default='1ST').
%        Its value can be one of the following:
%         - 1ST (first digit extraction)
%         - 2ND (second digit extraction)
%         - 3RD (third digit extraction)
%         - F2D (first-two digits extraction)
%         - L2D (last-two digits extraction)
% a    = A float [0.01,0.10] representing the statistical significance threshold for the Z-score tests (optional, default=0.05).
% ccf  = A boolean indicating whether to apply a continuity correction factor on Z-scores (optional, default=true).
% btv  = A boolean indicating whether to perform data transformation and validation (optional, default=true).
%        This input argument should be set to false only when data has already been transformed and validated by another function.
%
% [OUTPUT]
% tab  = A variable height table containing the extraction data, with the following columns:
%         - Digits: a vector of integers representing the ordered sequence of the significant digits.
%         - ThePL: a vector of floats representing the lower confidence intervals of the theoretical frequencies.
%         - TheP: a vector of floats representing the theoretical frequency of each significant digit.
%         - ThePU: a vector of floats representing the upper confidence intervals of the theoretical frequencies.
%         - TheF: a vector of floats representing the theoretical cumulative frequencies of the significant digits.
%         - EmpC: a vector of integers representing the number of occurrences of each significant digit in the sample.
%         - EmpP: a vector of floats representing the empirical frequency of each significant digit.
%         - EmpF: a vector of floats representing the empirical cumulative frequencies of the significant digits.
%         - Z: a vector of floats representing the Z-score of each significant digit.
%         - ZTest: a vector of booleans representing the Z-score test result of each significant digit.
% mad  = A 4-by-3 table containing the Mean Absolute Deviation test result, with the following columns:
%         - Conformity: a vector of strings representing the conformity level.
%         - Threshold: a vector of floats representing the conformity thresholds.
%         - MAD: a vector containing 3 NaN values and a float (the MAD value of the sample located at its respective conformity level).
% ssd  = A 4-by-3 table containing the Sum of Square Differences test result, with the following columns:
%         - Conformity: a vector of strings representing the conformity level.
%         - Threshold: a vector of floats representing the conformity thresholds.
%         - MAD: a vector containing 3 NaN values and a float (the SSD value of the sample located at its respective conformity level).

function [tab,mad,ssd] = benford_extract(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('dran','ALL',@(x)any(validatestring(x,{'ALL','NEG','POS'})));
        p.addOptional('ddec',2,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',10}));
        p.addOptional('extr','1ST',@(x)any(validatestring(x,{'1ST','2ND','3RD','F2D','L2D'})));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
        p.addOptional('ccf',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
        p.addOptional('btv',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    dran = res.dran;
    ddec = res.ddec;
    extr = res.extr;
    a = res.a;
    ccf = res.ccf;
    btv = res.btv;

    if (btv)
        data = benford_data(data,dran,ddec);
    end

    switch (nargout)
        case 1
            tab = benford_extract_internal(data,extr,a,ccf);

        case 2
            [tab,mad] = benford_extract_internal(data,extr,a,ccf);
            
        case 3
            [tab,mad,ssd] = benford_extract_internal(data,extr,a,ccf);

        otherwise
            error('Only up to 3 output arguments can be specified.');
    end

end

function [tab,mad,ssd] = benford_extract_internal(data,extr,a,ccf)

    switch (extr)
        case '1ST'
            [the_dgts,the_p,the_f,emp_c,emp_p,emp_f] = extract_first_digits(data,1);

        case '2ND'
            [the_dgts,the_p,the_f,emp_c,emp_p,emp_f] = extract_nth_digit(data,2);

        case '3RD'
            [the_dgts,the_p,the_f,emp_c,emp_p,emp_f] = extract_nth_digit(data,3);

        case 'F2D'
            [the_dgts,the_p,the_f,emp_c,emp_p,emp_f] = extract_first_digits(data,2);

        otherwise
            [the_dgts,the_p,the_f,emp_c,emp_p,emp_f] = extract_last_digits(data);
    end
    
    k = numel(the_dgts);
    n = numel(data);

    cv = norminv(1 - (a / 2));
    cv2 = cv ^ 2;

    diff_p = emp_p - the_p;
    diff_p_abs = abs(diff_p);
    v = (the_p .* (1 - the_p)) ./ n;
    
    if (ccf)
        wci_ccf = (4 .* the_p) + 2;
        wci_m = (2 .* n .* the_p) + cv2;
        wci_den = 2 * (n + cv2);
        wci_thr = cv2 - (1 / n) + (4 .* n .* the_p .* (1-the_p));

        the_p_lb = (wci_m - (cv .* sqrt(wci_thr + wci_ccf)) + 1) ./ wci_den;
        the_p_ub = (wci_m + (cv .* sqrt(wci_thr - wci_ccf)) + 1) ./ wci_den;

        z_ccf = 1 / (2 * n);
        z_num = diff_p_abs - min([diff_p_abs repmat(z_ccf,k,1)],[],2);
    else
        wci_den = 1 + (cv2 / n);
        wci_p = (the_p + (cv2 / (2 * n))) ./ wci_den;
        wci_moe = (cv / wci_den) .* sqrt(v + (cv2 / (4 * (n ^ 2))));

        the_p_lb = wci_p - wci_moe;
        the_p_ub = wci_p + wci_moe;

        z_num = diff_p_abs;
    end

    z = z_num ./ sqrt(v);
    z_test = z <= cv;

    tab = table(the_dgts,the_p_lb,the_p,the_p_ub,the_f,emp_c,emp_p,emp_f,z,z_test);
    tab.Properties.VariableNames = {'Digits' 'ThePL' 'TheP' 'ThePU' 'TheF' 'EmpC' 'EmpP' 'EmpF' 'Z' 'ZTest'};

    if (nargout >= 2)
        switch (extr)
            case '1ST'
                mad_thr = [0.006; 0.012; 0.015; Inf];

            case {'2ND' '3RD'}
                mad_thr = [0.008; 0.010; 0.012; Inf];

            otherwise
                mad_thr = [0.0012; 0.0018; 0.0022; Inf];
        end
        
        mad_val = mean(diff_p_abs);
        mad_idx = find(mad_thr >= mad_val,1);
        mad_res = [NaN((mad_idx - 1),1); mad_val; NaN((4 - mad_idx),1)];
        
        mad = table({'Close Conformity'; 'Acceptable Conformity'; 'Marginally Acceptable Conformity'; 'Nonconformity'},mad_thr,mad_res);
        mad.Properties.VariableNames = {'Conformity' 'Threshold' 'Value'};
    end
    
    if (nargout == 3)
        if (any(strcmp(extr,{'1ST' '2ND' '3RD'})))
            ssd_thr = [2; 25; 100; Inf];
        else
            ssd_thr = [2; 10; 50; Inf];
        end
        
        ssd_val = sum(diff_p .^ 2) * 10000;
        ssd_idx = find(ssd_thr > ssd_val,1);
        ssd_res = [NaN((ssd_idx - 1),1); ssd_val; NaN((4 - ssd_idx),1)];
        
        ssd = table({'Close Conformity'; 'Acceptable Conformity'; 'Marginally Acceptable Conformity'; 'Nonconformity'},ssd_thr,ssd_res);
        ssd.Properties.VariableNames = {'Conformity' 'Threshold' 'Value'};
    end

end

function [the_dgts,the_p,the_f,emp_c,emp_p,emp_f] = extract_first_digits(data,d)

    lo = 10 ^ (d - 1);
    hi = (10 ^ d) - 1;

    data = data(data >= lo);
    n = numel(data);

    the_dgts = (lo:hi).';
    the_p = log10(1 + (1 ./ the_dgts));
    the_f = [cumsum(the_p(1:end-1)); 1];

    emp_dgts = truncate((10 .^ ((floor(log10(data)) .* -1) + d - 1)) .* data);
    emp_c = histcounts(emp_dgts,[the_dgts; Inf]).';
    emp_p = emp_c ./ n;
    emp_f = [cumsum(emp_p(1:end-1)); 1];

end

function [the_dgts,the_p,the_f,emp_c,emp_p,emp_f] = extract_last_digits(data)

    data = data(data >= 10);
    n = numel(data);

    the_dgts = (0:99).';
    the_p = repmat(0.01,100,1);
    the_f = [cumsum(the_p(1:end-1)); 1];

    emp_dgts = mod(data,100);
    emp_c = histcounts(emp_dgts,[(0:99).'; Inf]).';
    emp_p = emp_c ./ n;
    emp_f = [cumsum(emp_p(1:end-1)); 1];

end

function [the_dgts,the_p,the_f,emp_c,emp_p,emp_f] = extract_nth_digit(data,d)

    data = data(data >= (10 ^ (d - 1)));
    n = numel(data);

    the_dgts = (0:9).';
    the_p = zeros(10,1);
    
    for i = 1:10
        k = ((10 ^ (d - 2)):((10 ^ (d - 1)) - 1)).';
        the_p(i) = sum(log10(1 + (1 ./ ((10 .* k) + the_dgts(i)))));
    end

    the_f = [cumsum(the_p(1:end-1)); 1];

    emp_k = max([floor(log10(data)) zeros(0,n,1)],[],2) + 2 - d;
    emp_dgts = floor(mod((data ./ (10 .^ (emp_k - 1))),10));
    emp_c = histcounts(emp_dgts,[the_dgts; Inf]).';
    emp_p = emp_c ./ n;
    emp_f = [cumsum(emp_p(1:end-1)); 1];

end

function x = truncate(x)

    x = str2double(cellstr(num2str(x)));
    x = x - rem(x,1);

end
