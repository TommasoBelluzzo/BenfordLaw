% [INPUT]
% data = A numeric array of n elements representing the sample on which the Number Duplication Analysis must be performed.
% ddec = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
%        No rounding is performed, the exceeding decimals are truncated as if they were not present.
%
% [OUTPUT]
% d10n = 
% d0n  = 
% d0p  = 
% d10p = 
% z    = 

function [d10n,d0n,d0p,d10p,z] = benford_duplication(varargin)

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

    if (nargout ~= 5)
        error('Only 5 output arguments can be specified.');
    end

    [d10n,d0n,d0p,d10p,z] = benford_duplication_internal(data,ddec);

end

function [d10n,d0n,d0p,d10p,z] = benford_duplication_internal(data,ddec)

    data = double(data(:));
    data = floor(data .* (10 ^ ddec)) ./ (10 ^ ddec);

    data_10n = data(data <= -10);

    if (numel(data_10n) > 1)
        [uni,~,uni_idx] = unique(data_10n);
        c = accumarray(uni_idx,1);

        d10n = table(uni,c);
        d10n.Properties.VariableNames = {'Value' 'Count'};
        
        d10n(d10n.Count < 2,:) = [];
        d10n = sortrows(d10n,[-2 1]);
    else
        d10n = table(zeros(0,1),zeros(0,1),'VariableNames',{'Value' 'Count'});
    end

    data_0n = data((data > -10) & (data < 0));
    
    if (numel(data_0n) > 1)
        [uni,~,uni_idx] = unique(data_0n);
        c = accumarray(uni_idx,1);

        d0n = table(uni,c);
        d0n.Properties.VariableNames = {'Value' 'Count'};
        
        d0n(d0n.Count < 2,:) = [];
        d0n = sortrows(d0n,[-2 1]);
    else
        d0n = table(zeros(0,1),zeros(0,1),'VariableNames',{'Value' 'Count'});
    end

    data_0p = data((data > 0) & (data < 10));
    
    if (numel(data_0p) > 1)
        [uni,~,uni_idx] = unique(data_0p);
        c = accumarray(uni_idx,1);

        d0p = table(uni,c);
        d0p.Properties.VariableNames = {'Value' 'Count'};
        
        d0p(d0p.Count < 2,:) = [];
        d0p = sortrows(d0p,[-2 -1]);
    else
        d0p = table(zeros(0,1),zeros(0,1),'VariableNames',{'Value' 'Count'});
    end

    data_10p = data(data >= 10);
    
    if (numel(data_10p) > 1)
        [uni,~,uni_idx] = unique(data_10p);
        c = accumarray(uni_idx,1);

        d10p = table(uni,c);
        d10p.Properties.VariableNames = {'Value' 'Count'};
        
        d10p(d10p.Count < 2,:) = [];
        d10p = sortrows(d10p,[-2 -1]);
    else
        d10p = table(zeros(0,1),zeros(0,1),'VariableNames',{'Value' 'Count'});
    end
    
    z = sum(data == 0);

end