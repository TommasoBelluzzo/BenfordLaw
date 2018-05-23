classdef BenfordData 
    %% Properties: Instance
    properties (GetAccess = public, SetAccess = private)
        OriginalSample;
        Sample;
        Table;
    end

    %% Constructor
    methods (Access = public)
        function this = BenfordData(os,s,t)
            this.OriginalSample = os;
            this.Sample = s;
            this.Table = t;
        end
    end
end1
