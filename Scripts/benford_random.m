% [INPUT]
% size = A vector of integers in which each element represents the size of a dimension of the output (optional, default=1);
% lim  = An integer, greater than or equal to 9, representing the maximum possible value to be generated (optional, default=9).
% prob = A string representing the method to use for distributing generated values over the whole range (optional, default='VAL').
%        Its value can be one of the following:
%         - MAG (each order of magnitude has the same probability)
%         - VAL (each value has the same probability)
%
% [OUTPUT]
% r    = An array of Benford's Law conforming random numbers.
%
% [NOTES]
% Credit goes to Jan Zdenek, the author of the original code.

function r = benford_random(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addOptional('size',1,@(x)validateattributes(x,{'numeric'},{'2d','nonempty','real','finite','integer','positive'}));
        p.addOptional('lim',9,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',9}));
        p.addOptional('prob','VAL',@(x)any(validatestring(x,{'MAG','VAL'})));
    end

    p.parse(varargin{:});

    res = p.Results;
    size = res.size;
    lim = res.lim;
    prob = res.prob;

    if (isscalar(size))
        size = repelem(size,1,2);
    end

    mag = floor(log10(lim) + 1);
    
    if ((mag == 1) && strcmp(prob,'MAG'))
        warning('Only one order of magnitude is present, probabilities will be handled on a per value basis.');
        prob = 'VAL';
    end
    
    r = benford_random_internal(size,lim,mag,prob);

end

function r = benford_random_internal(size,lim,mag,prob)

    x_len = prod(size);
    x_p = log10(1 + (1 ./ (1:9).'));
    x = datasample(1:9,x_len,'Weights',x_p);

    if strcmp(prob,'MAG')
        p = ones(mag,1) ./ mag;
    else
        p = [9; (9 .* (10 .^ ((2:mag-1).' - 1)))] ./ lim;
        p(end + 1) = 1 - sum(p);
    end

    idx = 1;
    stop = false(x_len,1);

    while (~all(stop))
        stop_neg = ~stop;
        x_max = (x .* 10) >= lim;
        stop(stop_neg & x_max) = true;
        
        stop_neg = ~stop;
        p_stop = rand(x_len,1) < p(idx);
        stop(stop_neg & p_stop) = true;

        x = arrayfun(@(n,s) add_digit(n,s,lim),x,stop);

        idx = idx + 1;
    end

    sizes = num2cell(size);
    r = reshape(x,sizes{:});

end

function x = add_digit(x,stop,lim)

    if (stop)
        return;
    end
    
    while (true)
        p_dgt = numel(num2str(x)) + 1;
        p_ran = (10 ^ (p_dgt - 2)):((10 ^ (p_dgt - 1)) - 1);
        p = arrayfun(@(d) sum(log10(1 + (1 ./ ((10 .* p_ran) + d)))),0:9).';

        out = (x * 10) + datasample(0:9,1,'Weights',p);

        if (out <= lim)
            x = out;
            break;
        end
    end

end