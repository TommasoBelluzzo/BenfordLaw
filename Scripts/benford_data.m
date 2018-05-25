% [INPUT]
% data = A numeric array representing the sample on which the Benford's Law analysis must be performed.
% d    = An integer [1,3] representing the number of first significant digits to analyse.
% ccf  = A boolean indicating whether to apply a continuity correction factor on Z-scores (optional, default=true).
%
% [OUTPUT]
% bd   = An instance of the BenfordData class that exposes the following properties:
%         - Digits: the number of first significant digits analysed.
%         - Sample: a vector containing the sample.
%         - SampleDigits: a vector containing the first significant digits of the sample.
%         - SampleOriginal: a vector containing the valid observations from which the sample has been produced.
%         - SecondOrder: a boolean indicating whether a second order analysis has been performed.
%         - Table: a kx8 table with the following columns:
%            > Digits: the ordered sequence of the first significant digits.
%            > Count: the number of occurrences of each first significant digit.
%            > Summation: the sum of the sample values grouped by their first significant digits, performed only for first order analyses.
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
        p.addOptional('ccf',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    d = res.d;
    ccf = res.ccf;

    data = double(data(:));
    data = data((data ~= 0) & isfinite(data));

    if (numel(unique(data)) < 10)
        error('The number of unique valid observations must be greater than or equal to 10.');
    end
    
    if (numel(data) < 250)
        warning('A minimum sample size of 250 unique observations is recommended in order to produce coherent analyses.');
    end
    
    bd = benford_data_internal(data,d,ccf);

end

function bd = benford_data_internal(data,d,ccf)

    the_dgts = ((10 ^ (d - 1)):((10 ^ d) - 1)).';
    the_dgts_len = numel(the_dgts);
    the_p = log10(1 + (1 ./ the_dgts));
    the_f = [cumsum(the_p(1:end-1)); 1];
    the_su_p = repmat((1 / the_dgts_len),the_dgts_len,1);
    the_su_f = [cumsum(the_su_p(1:end-1)); 1];

    fo = abs(data);
    fo_dgts = truncate((10 .^ ((floor(log10(fo)) .* -1) + d - 1)) .* fo);
    [fo_a,fo_p,fo_f,fo_z] = analyse_standard(fo_dgts,the_dgts,the_p,ccf);
    
    fo_tab = table(the_dgts,fo_a,the_p,fo_p,the_f,fo_f,fo_z);
    fo_tab.Properties.VariableNames = {'Digits' 'Amount' 'TheP' 'EmpP' 'TheF' 'EmpF' 'Z'};
    
    so = generate_second_order(fo);
    so_dgts = truncate((10 .^ ((floor(log10(so)) .* -1) + d - 1)) .* so);
    [so_a,so_p,so_f,so_z] = analyse_standard(so_dgts,the_dgts,the_p,ccf);

    so_tab = table(the_dgts,so_a,the_p,so_p,the_f,so_f,so_z);
    so_tab.Properties.VariableNames = {'Digits' 'Amount' 'TheP' 'EmpP' 'TheF' 'EmpF' 'Z'};

    [su_a,su_p,su_f,su_z] = analyse_summation(the_dgts,the_p,fo,fo_dgts,ccf);

    su_tab = table(the_dgts,su_a,the_su_p,su_p,the_su_f,su_f,su_z);
    su_tab.Properties.VariableNames = {'Digits' 'Amount' 'TheP' 'EmpP' 'TheF' 'EmpF' 'Z'};
    
    bd = BenfordData(d,data,fo,fo_dgts,fo_tab,so,so_dgts,so_tab,su_tab);

end

function [a,p,f,z] = analyse_standard(emp_dgts,the_dgts,the_p,ccf)

    n = numel(emp_dgts);

    a = histcounts(emp_dgts,[the_dgts; Inf]).';
    p = a ./ n;
    f = [cumsum(p(1:end-1)); 1];
    z = calculate_z(n,p,the_p,ccf);

end

function [s,p,f,z] = analyse_summation(the_dgts,the_p,fo,fo_dgts,ccf)

    [dgts_uni,~,dgts_uni_idx] = unique(fo_dgts); 
    su_tab = [dgts_uni accumarray(dgts_uni_idx,fo)];
    su_tab_idx = ismember(the_dgts,dgts_uni);

    s = zeros(numel(the_dgts),1);
    s(su_tab_idx) = su_tab(:,2);
    p = s ./ sum(s);
    f = [cumsum(p(1:end-1)); 1];
    z = calculate_z(numel(fo_dgts),p,the_p,ccf);

end

function z = calculate_z(n,emp_p,the_p,ccf)

    z_num = abs(emp_p - the_p);
    z_den = sqrt((the_p .* (1 - the_p)) ./ n);
    
    if (ccf)
        ccf_val = 1 / (2 * n);

        ccf_idx = z_num > ccf_val;
        z_num(ccf_idx) = z_num(ccf_idx) - ccf_val;
    end
    
    z = z_num ./ z_den;

end

function data = generate_second_order(data)

    data = sort(data);
    data = abs(round(data(2:end) - data(1:end-1),10));
    data = data(data ~= 0);

end

function x = truncate(x)

    x = x - rem(x,1);

end