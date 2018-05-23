% [INPUT]
% data = A numeric array representing the sample on which the Benford's Law analysis must be performed.
% d    = An integer [1,3] representing the number of first significant digits to analyse.
% sec  = A boolean indicating whether to perform a second order analysis (optional, default=false).
% ccf  = A boolean indicating whether to apply a continuity correction factor on Z-scores (optional, default=true).
%
% [OUTPUT]
% bd   = An instance of the BenfordData class that exposes the following properties:
%         - OriginalSample: a vector containing the valid observations.
%         - Sample: a vector containing the post-processed valid observations.
%         - Table: a kx7 table with the following columns:
%            > Digits: the ordered sequence of the first significant digits.
%            > Count: the number of occurrences of each first significant digit.
%            > TheP: the theoretical frequency of each first significant digit.
%            > EmpP: the empirical frequency of each first significant digit.
%            > TheF: the theoretical cumulative frequency of the first significant digits.
%            > EmpF: the empirical cumulative frequency of the first significant digits.
%            > Z: the Z-score of each first significant digit.

function bd = benford_data(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addRequired('d',@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',1,'<=',3}));
        p.addOptional('sec',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
        p.addOptional('ccf',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    d = res.d;
    sec = res.sec;
    ccf = res.ccf;

    data = double(data(:));
    data = data((data ~= 0) & isfinite(data));
    
    data_len = numel(data);
    
    if (data_len < 10)
        error('The number of valid observations must be greater than or equal to 10.');
    end
    
    if (data_len < 250)
        warning('A minimum sample size of 250 unique observations is recommended in order to produce coherent analyses.');
    end
    
    bd = benford_data_internal(data,d,sec,ccf);

end

function bd = benford_data_internal(data,d,sec,ccf)

    data_orig = data;

    data = abs(data);

    if (sec)
        data = sort(data);
        data = abs(round(data(2:end) - data(1:end-1),10));
        data = data(data ~= 0);
    end
    
    data_len = numel(data);

    the_dgts = ((10 ^ (d - 1)):((10 ^ d) - 1)).';
    the_p = log10(1 + (1 ./ the_dgts));
    the_f = cumsum(the_p);
    
    emp_dgts = (10 .^ ((floor(log10(data)) .* -1) + d - 1)) .* data;
    emp_dgts = emp_dgts - rem(emp_dgts,1);
    
    emp_hist = histcounts(emp_dgts,[the_dgts; Inf]).';
    emp_p = emp_hist ./ data_len;
    emp_f = cumsum(emp_p);

    z_num = abs(emp_p - the_p);
    z_den = sqrt((the_p .* (1 - the_p)) ./ data_len);
    
    if (ccf)
        ccf_val = 1 / (2 * data_len);
        ccf_idx = z_num > ccf_val;

        z_num(ccf_idx) = z_num(ccf_idx) - ccf_val;
    end
    
    emp_z = z_num ./ z_den;

    tab = table(the_dgts,emp_hist,the_p,emp_p,the_f,emp_f,emp_z);
    tab.Properties.VariableNames = {'Digits' 'Count' 'TheP' 'EmpP' 'TheF' 'EmpF' 'Z'};
    
    bd = BenfordData(data_orig,data,tab);

end