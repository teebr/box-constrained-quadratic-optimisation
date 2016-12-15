classdef qpBoxConstTest < qpBox
    %An example for a constant QP problem
    
    %H and f are both constant in this case, so just set them as properties
    %and have the hCalc and fCalc functions call the property
    properties
        H
        f
    end
    
    methods
        function obj = qpBoxConstTest(varargin)
            setProperties(obj,nargin,varargin{:});
        end
        
        function H = hCalc(obj,varargin) %need varargin even if no args
            H = obj.H;
        end
        
        function f = fCalc(obj,varargin) %need varargin even if no args
            f = obj.f;
        end
    end
    
    methods (Access=protected)
        function num = getNumInputsImpl(~)
            num = 0;
        end
    end
end