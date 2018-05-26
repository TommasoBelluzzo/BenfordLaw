classdef BenfordData 
    %% Properties: Instance
    properties (GetAccess = public, SetAccess = private)
        Digits;
        Sample;
        FirstOrderData;
        FirstOrderDigits;
        FirstOrderTable;
        SecondOrderData;
        SecondOrderDigits;
        SecondOrderTable;
        Mantissae;
        Summation;
    end

    %% Constructor
    methods (Access = public)
        function this = BenfordData(d,data,fo,fo_dgts,fo_tab,so,so_dgts,so_tab,mant,su)
            this.Digits = d;
            this.Sample = data;
            this.FirstOrderData = fo;
            this.FirstOrderDigits = fo_dgts;
            this.FirstOrderTable = fo_tab;
            this.SecondOrderData = so;
            this.SecondOrderDigits = so_dgts;
            this.SecondOrderTable = so_tab;
            this.Mantissae = mant;
            this.Summation = su;
        end
    end
end

