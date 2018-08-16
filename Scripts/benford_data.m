% [INPUT]
% data = A numeric array representing the sample to transform and validate.
% dran = A string representing the range of values to consider (optional, default='ALL').
%        Its value can be one of the following:
%         - ALL (all values)
%         - NEG (only negative values)
%         - POS (only positive values)
% ddec = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
%        No rounding is performed, the exceeding decimals are truncated as if they were not present.
%
% [OUTPUT]
% data = A numeric vector representing the transformed and validated sample.

function data = benford_data(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('dran','ALL',@(x)any(validatestring(x,{'ALL','NEG','POS'})));
        p.addOptional('ddec',2,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',10}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    dran = res.dran;
    ddec = res.ddec;

    data = benford_data_internal(data,dran,ddec);

end

function data = benford_data_internal(data,dran,ddec)

    data = double(data(:));
    data = floor(data .* (10 ^ ddec));
    
    switch (dran)
        case 'NEG'
            data = abs(data(isfinite(data) & (data < 0)));
    
        case 'POS'
            data = data(isfinite(data) & (data > 0));
            
        otherwise
            data = abs(data(isfinite(data) & (data ~= 0)));
    end

    if (numel(unique(data)) < 50)
        error('The number of unique valid observations in the sample must be greater than or equal to 50.');
    end
    
    if (numel(data) < 1000)
        warning('A minimum sample size of 1000 valid observations is recommended in order to produce a consistent analysis.');
    end

end
