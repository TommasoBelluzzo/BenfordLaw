% [INPUT]
% data = A numeric array representing the sample on which the Benford's Law analysis must be performed.

function benford_analyse(varargin)

    persistent p;

    if (isempty(p))
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    
    benford_analyse_internal(data);

end

function r = benford_analyse_internal(size,lim,mag,prob)

    x_len = prod(size);
    x = datasample(1:9,x_len,'Weights',log10(1 + (1 ./ (1:9).'))).';

    if strcmp(prob,'MAG')
        p_sam = ones(mag,1) ./ mag;
    else
        p_sam = [9; (9 .* (10 .^ ((2:mag-1).' - 1)))] ./ lim;
        p_sam(end+1) = 1 - sum(p_sam);
    end
    
    p_tab_len = numel(num2str(lim)) - 1;
    p_tab = NaN(10,p_tab_len);
    
    for i = 1:p_tab_len
        p_dgt = i + 1;
        p_ran = (10 ^ (p_dgt - 2)):((10 ^ (p_dgt - 1)) - 1);
        p_tab(:,i) = arrayfun(@(d) sum(log10(1 + (1 ./ ((10 .* p_ran) + d)))),0:9).';
    end
    
    k = 1;
    stop = false(x_len,1);

    while (~all(stop))
        stop(~stop & (((x .* 10) >= lim) | (rand(x_len,1) < p_sam(k)))) = true;
        x(~stop) = arrayfun(@(n) add_digit(n,lim,p_tab),x(~stop));

        k = k + 1;
    end

    sizes = num2cell(size);
    r = reshape(x,sizes{:});

end

function x = add_digit(x,lim,p_tab)

    while (true)
        d = numel(num2str(x)) + 1;
        p = p_tab(:,(d - 1));

        out = (x * 10) + datasample(0:9,1,'Weights',p);

        if (out <= lim)
            x = out;
            break;
        end
    end

end