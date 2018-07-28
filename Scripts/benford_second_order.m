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
% tab  = A 90-by-6 table containing the second order data, with the following columns:
%         - Digits: a vector of integers representing the ordered sequence of the first-two digits.
%         - TheP: a vector of floats representing the theoretical frequencies of the first-two digits.
%         - TheF: a vector of floats representing the theoretical cumulative frequencies of the first-two digits.
%         - EmpC: a vector of integers representing the number of occurrences of each pair of digits.
%         - EmpP: a vector of floats representing the empirical frequency of each pair of digits.
%         - EmpF: a vector of floats representing the empirical cumulative frequencies of each pair of digits.

function tab = benford_second_order(varargin)

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

	tab = benford_second_order_internal(data);

end

function tab = benford_second_order_internal(data)

	data = sort(data);
    data = diff(data);
    data = data(data ~= 0);
    
    n = numel(data);

    the_dgts = (10:99).';
    the_p = log10(1 + (1 ./ the_dgts));
    the_f = [cumsum(the_p(1:end-1)); 1];

    emp_dgts = (10 .^ ((floor(log10(data)) .* -1) + 1)) .* data;
    emp_dgts = str2double(cellstr(num2str(emp_dgts)));
    emp_dgts = emp_dgts - rem(emp_dgts,1);

    emp_c = histcounts(emp_dgts,[the_dgts; Inf]).';
    emp_p = emp_c ./ n;
    emp_f = [cumsum(emp_p(1:end-1)); 1];

    tab = table(the_dgts,the_p,the_f,emp_c,emp_p,emp_f);
    tab.Properties.VariableNames = {'Digits' 'TheP' 'TheF' 'EmpC' 'EmpP' 'EmpF'};

end