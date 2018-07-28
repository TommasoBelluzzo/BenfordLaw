% [INPUT]
% data = A numeric array of n elements representing the original (non post-processed) sample on which the Summation Analysis must be performed.
% ddec = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
%        No rounding is performed, the exceeding decimals are truncated as if they were not present.
%
% [OUTPUT]
% tab  = A 90-by-7 table containing the summation data, with the following columns:
%         - Digits: a vector of integers representing the ordered sequence of the first-two digits.
%         - Summation: a vector of floats representing the summation of the sample elements having the same first-two digits.
%         - AES: a vector of floats representing the absolute excess summations.
%         - TheP: a vector of floats representing the theoretical frequencies of the first-two digits.
%         - TheF: a vector of floats representing the theoretical cumulative frequencies of the first-two digits.
%         - EmpP: a vector of floats representing the empirical frequency of each pair of digits.
%         - EmpF: a vector of floats representing the empirical cumulative frequencies of each pair of digits.

function tab = benford_summation(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('ddec',2,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',10}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    ddec = res.ddec;

    tab = benford_summation_internal(data,ddec);

end

function tab = benford_summation_internal(data,ddec)

    exp = 10 ^ ddec;

    data = double(data(:));
    data = data(isfinite(data) & (data >= 10));
    data = floor(data .* exp) ./ exp;

    the_dgts = (10:99).';
    the_p = repmat(1/90,90,1);
    the_f = [cumsum(the_p(1:end-1)); 1];

    emp_dgts = (10 .^ ((floor(log10(data)) .* -1) + 1)) .* data;
    emp_dgts = str2double(cellstr(num2str(emp_dgts)));
    emp_dgts = emp_dgts - rem(emp_dgts,1);

    [emp_dgts_uni,~,emp_dgts_idx] = unique(emp_dgts); 
    s_tab = [emp_dgts_uni accumarray(emp_dgts_idx,data)];

    s = zeros(90,1);
    s(ismember(the_dgts,emp_dgts_uni)) = s_tab(:,2);
    
    aes = abs(s - round(mean(s),2));

    emp_p = s ./ sum(s);
    emp_f = [cumsum(emp_p(1:end-1)); 1];

    tab = table(the_dgts,s,aes,the_p,the_f,emp_p,emp_f);
    tab.Properties.VariableNames = {'Digits' 'Summation' 'AES' 'TheP' 'TheF' 'EmpP' 'EmpF'};

end