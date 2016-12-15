classdef qpBoxFcnTest < qpBox
    %An example for a constant QP problem
    
    properties(Nontunable)
        N = 10;
    end
    
    methods
        function obj = qpBoxFcnTest(varargin)
            setProperties(obj,nargin,varargin{:});
        end
        
        function H = hCalc(obj,a,~)
            H = a * eye(obj.N) -diag(ones(obj.N-1,1),1) -diag(ones(obj.N-1,1),-1);
        end
        
        function f = fCalc(obj,a,b)
            f = b*sin(a*(1:obj.N)');
        end
    end
    
    methods (Access=protected)
        function num = getNumInputsImpl(~)
            num = 2;
        end
        
        function validatePropertiesImpl(obj)
            superValidate(obj) %must be called to set up qpBox
            %now add things specific to this instance
            
            if obj.N ~= 10
                disp('warning: demo won''t work')
            end
        end
    end
end