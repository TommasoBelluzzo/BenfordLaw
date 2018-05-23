% [INPUT]
% size = A vector of integers in which each element represents the size of a dimension of the output (optional, default=1);
% lim  = An integer, greater than or equal to 9, representing the maximum possible value to be generated (optional, default=9).
% rgen = A string representing the method to use for random numbers generation (optional, default='NAT').
%        Its value can be one of the following:
%         - ART (artificial, random numbers are used to increase the count of the distribution bins according to the Benford's Law frequencies)
%         - NAT (natural, random numbers between 0 and 1 are used to create a log-uniform distribution)
% prob = A string representing the method to use for distributing generated values over the whole range (optional, default='VAL').
%        It should only be specified if the artificial random numbers generation is being used. Its value can be one of the following:
%         - MAG (each order of magnitude has the same probability)
%         - VAL (each value has the same probability)
%
% [OUTPUT]
% x    = An array of Benford's Law conforming random numbers.
%
% [NOTES]
% Credit goes to Jan Zdenek, the author of the original code.

function r = benford_random(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addOptional('size',1,@(x)validateattributes(x,{'numeric'},{'2d','nonempty','real','finite','integer','positive'}));
        p.addOptional('lim',9,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',9}));
        p.addOptional('rgen','NAT',@(x)any(validatestring(x,{'ART','NAT'})));
        p.addOptional('prob','VAL',@(x)any(validatestring(x,{'MAG','VAL'})));
    end

    p.parse(varargin{:});

    res = p.Results;
    size = res.size;
    lim = res.lim;
    rgen = res.rgen;
    prob = res.prob;

    if ((nargin == 4) && strcmp(rgen,'NAT'))
        error('The ''prob'' parameter should only be specified if the artificial random numbers generation is being used.');
    end

    if (isscalar(size))
        size = repelem(size,1,2);
    end

    mag = floor(log10(lim) + 1);
    
    if ((mag == 1) && strcmp(prob,'MAG'))
        warning('Only one order of magnitude is present, probabilities will be handled on a per value basis.');
        prob = 'VAL';
    end
    
    r = benford_random_internal(size,lim,prob,rgen);

end

function r = benford_random_internal(size,lim,prob,rgen)

    if (strcmp(rgen,'NAT'))
        x = floor(exp(log(lim) * rand(size)));
    else
        stoppers = zeros(1,prod(size));
        prob_stops = getprobs(floor(log10(lim) + 1), 1, lim,prob);
        i = 1;

        x = arrayfun(@(x) setfirstdigit(), stoppers);

        while prod(stoppers) == 0
            stoppers = arrayfun(@(n, s) shallstop(n, s, 1, lim, prob_stops(i)), x, stoppers);
            x = arrayfun(@(s, n) addnextdigit(n, s, lim), stoppers, x);
            i = i+1;
        end
    end

    sizes = num2cell(size);
    r = reshape(x,sizes{:});

end

function digit = setfirstdigit()
distribution = [0.3010 0.4771 0.6020 0.6980 0.7772 0.8441 0.9021 0.9533 1.0];
over_rand_threshold = find(rand(1) < distribution);
digit = over_rand_threshold(1);
end

function number = addnextdigit(number, stopper, max)
% Adds another digit to the number unless the generation of the number has
% been terminated.
if stopper == 1
    return;
end
addnext = 1;
while addnext == 1
    digit = floor(10*rand(1));
    if number*10 + digit <= max % Makes sure that the new number doesn't exceed max value
        number = number*10 + digit;
        addnext = 0;
    end
end
end

function stopper = shallstop(number, stopper, min, max, prob_stop)
% Determines whether the generation of the number will be stopped or not.
% The number must exceed min value and must not exceed max value, and its
% generation can be terminated with a certain probability.
if stopper == 1
    return;
end
if number < min
    stopper = 0;
elseif number * 10 >= max
    stopper = 1;
else
    if rand(1) < prob_stop
        stopper = 1;
    else
        stopper = 0;
    end
end
end

function probabilities = getprobs(num_of_digits, min, max, prob_method)
% Creates an array of probabilities of termination of generation process
% for each order of magnitude.
probabilities = zeros(1, num_of_digits);
min_num_of_digits = floor(log10(min) + 1);

for i = min_num_of_digits:num_of_digits
    if strcmp(prob_method,'MAG')
        probabilities(i) =  (1 / (num_of_digits - min_num_of_digits + 1));
    else
        if i == min_num_of_digits
            probabilities(i) = (10^(i) - min) / (max - min + 1);
        elseif i == num_of_digits
            probabilities(i) = ((9 * 10^(i-1)) - (10^i - max)) / (max - min + 1);
        else
            probabilities(i) =  (9 * 10^(i-1)) / (max - min + 1);
        end
    end
end
end