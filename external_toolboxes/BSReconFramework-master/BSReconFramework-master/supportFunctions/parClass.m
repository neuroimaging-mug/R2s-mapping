classdef parClass < handle
    properties
        par = [];
    end
    methods
        function obj = parClass(par)
            obj.par = par;
        end
    end
end