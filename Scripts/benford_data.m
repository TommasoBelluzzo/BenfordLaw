% [INPUT]
% data = A numeric array representing the sample to transform and validate.
% ran  = A string representing the range of values to consider (optional, default='ALL').
%        Its value can be one of the following:
%         - ALL (all values)
%         - NEG (only negative values)
%         - POS (only positive values)
% dec  = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
%
% [OUTPUT]
% data = A numeric vector representing the transformed and validated sample.

function data = benford_data(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('ran','ALL',@(x)any(validatestring(x,{'ALL','NEG','POS'})));
        p.addOptional('dec',2,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',10}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    ran = res.ran;
    dec = res.dec;

    data = benford_data_internal(data,ran,dec);

end

function data = benford_data_internal(data,ran,dec)

    data = double(data(:));
    data = round(data .* (10 ^ dec),0);
    
    switch (ran)
        case 'NEG'
            data = abs(data(isfinite(data) & (data < 0)));
    
        case 'POS'
            data = data(isfinite(data) & (data > 0));
            
        otherwise
            data = abs(data(isfinite(data) & (data ~= 0)));
    end

    if (numel(unique(data)) < 10)
        error('The number of unique valid observations in the sample must be greater than or equal to 50.');
    end
    
    if (numel(data) < 1000)
        warning('A minimum sample size of 1000 valid observations is recommended in order to produce a coherent analysis.');
    end

end

% [INPUT]
% data = A numeric array representing the sample on which the Benford's Law Analysis must be performed.
% d    = An integer [1,3] representing the number of first significant digits to analyse.
% ccf  = A boolean indicating whether to apply a continuity correction factor on Z-scores (optional, default=true).
%
% [OUTPUT]
% bd   = An instance of the BenfordData class that exposes the following properties:
%         - Digits = The number of first significant digits analysed.
%         - Sample = A vector containing the original valid observations.
%         - FirstOrderData = A vector containing the first order sample.
%         - FirstOrderDigits = A vector containing the first significant digits of the first order sample.
%         - FirstOrderTable = A n-by-7 table containing the first order test data. It has the following columns:
%            > Digits: the ordered sequence of the first significant digits.
%            > Amount: the number of occurrences of each first significant digit.
%            > TheP: the theoretical frequency of each first significant digit.
%            > EmpP: the empirical frequency of each first significant digit.
%            > TheF: the theoretical cumulative frequency of the first significant digits.
%            > EmpF: the empirical cumulative frequency of the first significant digits.
%            > Z: the Z-score of each first significant digit.
%         - SecondOrderData = A vector containing the second order sample.
%         - SecondOrderDigits = A vector containing the first significant digits of the second order sample.
%         - SecondOrderTable = A n-by-7 table containing the second order test data. It has the following columns:
%            > Digits: the ordered sequence of the first significant digits.
%            > Amount: the number of occurrences of each first significant digit.
%            > TheP: the theoretical frequency of each first significant digit.
%            > EmpP: the empirical frequency of each first significant digit.
%            > TheF: the theoretical cumulative frequency of the first significant digits.
%            > EmpF: the empirical cumulative frequency of the first significant digits.
%            > Z: the Z-score of each first significant digit.
%         - Summation = A n-by-7 table containing the summation test data. It has the following columns:
%            > Digits: the ordered sequence of the first significant digits.
%            > Amount: the sum of the sample values grouped by their respective first significant digits.
%            > TheP: the theoretical frequency of each sum.
%            > EmpP: the empirical frequency of each sum.
%            > TheF: the theoretical cumulative frequencies of the sums.
%            > EmpF: the empirical cumulative frequencies of the sums.
%            > Z: the Z-score of each sum.
% 
% function bd = benford_data(varargin)
% 
%     persistent p;
% 
%     if (isempty(p))
%         p = inputParser();
%         p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
%         p.addRequired('d',@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',1,'<=',3}));
%         p.addOptional('ccf',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
%     end
% 
%     p.parse(varargin{:});
% 
%     res = p.Results;
%     data = res.data;
%     d = res.d;
%     ccf = res.ccf;
% 
%     data = double(data(:));
%     data = data((data ~= 0) & isfinite(data));
% 
%     if (numel(unique(data)) < 10)
%         error('The number of unique valid observations must be greater than or equal to 10.');
%     end
%     
%     if (numel(data) < 250)
%         warning('A minimum sample size of 250 unique observations is recommended in order to produce coherent analyses.');
%     end
%     
%     bd = benford_data_internal(data,d,ccf);
% 
% end
% 
% function bd = benford_data_internal(data,d,ccf)
% 
%     the_dgts = ((10 ^ (d - 1)):((10 ^ d) - 1)).';
%     the_dgts_len = numel(the_dgts);
%     the_p = log10(1 + (1 ./ the_dgts));
%     the_f = [cumsum(the_p(1:end-1)); 1];
%     the_su_p = repmat((1 / the_dgts_len),the_dgts_len,1);
%     the_su_f = [cumsum(the_su_p(1:end-1)); 1];
% 
%     fo = abs(data);
%     fo_dgts = truncate((10 .^ ((floor(log10(fo)) .* -1) + d - 1)) .* fo);
%     [fo_a,fo_p,fo_f,fo_z] = analyse_order(fo_dgts,the_dgts,the_p,ccf);
%     fo_tab = table(the_dgts,fo_a,the_p,fo_p,the_f,fo_f,fo_z);
%     fo_tab.Properties.VariableNames = {'Digits' 'Amount' 'TheP' 'EmpP' 'TheF' 'EmpF' 'Z'};
%     
%     so = generate_second_order(data);
%     so_dgts = truncate((10 .^ ((floor(log10(so)) .* -1) + d - 1)) .* so);
%     [so_a,so_p,so_f,so_z] = analyse_order(so_dgts,the_dgts,the_p,ccf);
%     so_tab = table(the_dgts,so_a,the_p,so_p,the_f,so_f,so_z);
%     so_tab.Properties.VariableNames = {'Digits' 'Amount' 'TheP' 'EmpP' 'TheF' 'EmpF' 'Z'};
%     
%     [su_a,su_p,su_f,su_z] = analyse_summation(the_dgts,the_p,fo,fo_dgts,ccf);
%     su = table(the_dgts,su_a,the_su_p,su_p,the_su_f,su_f,su_z);
%     su.Properties.VariableNames = {'Digits' 'Amount' 'TheP' 'EmpP' 'TheF' 'EmpF' 'Z'};
%     
%     bd = BenfordData(d,data,fo,fo_dgts,fo_tab,so,so_dgts,so_tab,mant,su);
% 
% end
% 
% function [s,p,f,z] = analyse_summation(the_dgts,the_p,fo,fo_dgts,ccf)
% 
%     k = numel(the_dgts);
%     n = numel(fo_dgts);
% 
%     [dgts_uni,~,dgts_uni_idx] = unique(fo_dgts); 
%     su_tab = [dgts_uni accumarray(dgts_uni_idx,fo)];
%     su_tab_idx = ismember(the_dgts,dgts_uni);
% 
%     s = zeros(k,1);
%     s(su_tab_idx) = su_tab(:,2);
%     p = s ./ sum(s);
%     f = [cumsum(p(1:end-1)); 1];
%     z = calculate_z(n,p,the_p,ccf);
% 
% end
% 
% function data = generate_second_order(data)
% 
%     data = diff(sort(data));
%     data = data(data ~= 0);
%     data = abs(data);
% 
% end
% 
