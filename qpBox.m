classdef (Abstract) qpBox < matlab.System
    % Solve a box-constrained quadratic minimisation problem. This is an
    % abstract class: it needs to be inherited to add functions to compute
    % the H matrix and f vector. These must be called hCalc and fCalc
    % respectively, and can be a function of the inputs and any of
    % properties you define. Additionally, the getNumInputsImpl method must
    % be defined (this is normally defined by the stepImpl() function but
    % here that is left as varargin so it will work with any number of
    % inputs)
    %
    % The upper and lower bounds and f vector must all be column vectors.
    % Some basic checks are performed for dimensions.
    % There are two outputs: the solution vector, and a status value. The
    % status value takes the following values:
    %   0: failed (quadratic function value is increasing)
    %   1: converging but did not reach function tolerance: increase
    %      tolerance or maximum number of iterations
    %   2: Optimisation solved successfully
    %
    % This was originally written to as part of a controller running as
    % generated code from Simulink. If you have a specific problem, it
    % might be more efficient to copy all the solver methods and write the
    % H and f functions int to the code directly.
    properties
        tolFun = 1e-8; %stopping tolerance
        gamma = 0.8; %step size damping
        maxIter = 25;
        lb = 0;
        ub = 1;
    end
    
    properties (Access = protected)
        gammaStop %this is updated in setup (basically a dependent)
        x0 = []; %this is updated in setup (basically a dependent)
        numElements = 1; %number of elements
    end
    
    methods (Access = protected)
        function validatePropertiesImpl(obj)
            %call core validation function. This function lets you easily
            %add your own validation checks on top of these. Just overwrite
            %the validateProperties function in your class and all
            %superValidate before adding whatever else
            superValidate(obj)
        end
        
        function setupImpl(obj,varargin)
            %see validateProperties- same applies.
            superSetup(obj,varargin{:})
        end
        
        function [u,statusFlag] = stepImpl(obj,varargin)
            %get H and f values, then solve problem
            H = obj.hCalc(varargin{:});
            if any(eig(H) < 0)
                disp('warning: H should be positive definite') %can't use 
            end
            f = obj.fCalc(varargin{:});
            [u,statusFlag] = obj.mainSolver(H,f);
        end %StepImpl
        
        function num = getNumOutputsImpl(~)
            % optimisation vector and status flag
            num = 2;
        end
        
        function [name1, name2] = getOutputNamesImpl(~)
            name1 = 'u';
            name2 = 'statusFlag';
        end
        
        function [fz1, fz2] = isOutputFixedSizeImpl(~)
            %Both outputs are always fixed-size
            fz1 = true;
            fz2 = true;
        end
        
        function [sz1 , sz2] = getOutputSizeImpl(obj)
            sz1 = [obj.numElements,1];
            sz2= [1 1]; %scalar value for errorFlag
        end
        
        function [dt1, dt2] = getOutputDataTypeImpl(~)
            dt1 = 'double';
            dt2 = 'double';
        end
        
        function [cp1, cp2] = isOutputComplexImpl(~)
            cp1 = false;
            cp2 = false;
        end
        
        function superValidate(obj)
            %setup for this function
            if ~iscolumn(obj.lb)
                error('Lower bound must be a column vector')
            end
            if ~iscolumn(obj.ub)
                error('Upper bound must be a column vector')
            end
            if numel(obj.lb) ~= numel(obj.ub)
                error('Upper and lower bounds must be the same length')
            end
            if any(obj.lb >= obj.ub)
                error('bounds are inconsistent: all elements of lower bound must be less than upper bound')
            end
            if ~isscalar(obj.tolFun) || obj.tolFun < 0
                error('tolFun must be a positive scalar value')
            end
            if ~isscalar(obj.maxIter) || obj.maxIter <= 0
                error('maxIter must be a positive scalar value')
            end
            if ~isscalar(obj.gamma) || obj.gamma <= 0 || obj.gamma >= 1
                error('gamma must be a scalar between 0 and 1')
            end
        end
        
        function superSetup(obj,varargin)
            obj.gammaStop = obj.gamma^5; %allows 5 attempts
            obj.numElements = numel(obj.lb);
            
            %determine initial guess for x0
            x = 0.5*(obj.lb + obj.ub);
            uInf = obj.ub == inf;
            lInf = obj.lb == -inf;
            x(uInf) = obj.lb(uInf)*1.2;
            x(lInf) = obj.ub(lInf)*0.8;
            x(lInf & uInf) = 0;

            x = obj.moveInBounds(x,obj.lb,obj.ub); %make sure inside of bounds
            obj.x0 = x;
            
            %check H and f dimensions
            H = obj.hCalc(varargin{:});
            f = obj.fCalc(varargin{:});
            if size(H,1) ~= obj.numElements || size(H,2) ~= obj.numElements
                error('H must be a square matrix with a length equal to f and the bound vectors')
            end
            if ~iscolumn(f) || numel(f) ~= obj.numElements
                error('f must be a column vector with a length equal to H and the bound vectors')
            end
        end
        
        function [x,statusFlag] = mainSolver(obj,H,f)
            x = obj.x0;
            statusFlag = 1; %default: not failed, but not reached requested tolerance
            gam = obj.gamma;
            fVal = obj.fQuad(H,f,x); %initial function value.
            
            %Aim to solve D^-2(x)g(x) = 0 using PCG method.
            for n = 1:obj.maxIter
                %calculate g and D
                grad = 2*H*x + f; %gradient of quadratic function
                D2 = obj.scalingMatrix(x,grad,obj.lb,obj.ub);
                
                %solve new step (eqns 9.9 - to 9.11, but simplified)
                GJv = diag(abs(grad)); %diag(g) * Jv (Jv is jacobian of g, i.e. diag*sign(g))
                A = H + GJv*D2; %D*GJv*D but D is diagonal and can direrctly calculate D*D
                s = A\-grad;
                
                %calculate next step:
                [xHat,dfval,fValHat] = obj.nextStep(x,s,gam,obj.lb,obj.ub,fVal,H,f);
                
                %check new step is okay, adjust step size if necessary
                while dfval > 0 && gam > obj.gammaStop %dfval > 0: function value is increasing
                    gam = gam*obj.gamma; %update damping
                    [xHat,dfval,fValHat] = obj.nextStep(x,s,gam,obj.lb,obj.ub,fVal,H,f);
                end
                
                if -dfval < obj.tolFun
                    %finished
                    statusFlag = 2;
                    break
                    
                elseif gam <= obj.gammaStop
                    % give up reducing step size - function doesn't appear to be decreasing
                    statusFlag = 0;
                    break
                else
                    %new step is OK, so update values
                    x = xHat;
                    fVal = fValHat;
                    gam = obj.gamma; %reset damping
                end
            end
            %complete
        end
        
    end
    
    methods (Static)
        % SUB-ROUTINES FOR OPTIMISATION
        
        function q = fQuad(H,f,x)
            %evaluate quadratic function. Only one line but anonymous function is slow.
            q = x'*H*x + f'*x;
        end
        
        function D2 = scalingMatrix(x,g,lb,ub)
            %D scaling matrix
            
            v = zeros(size(x));
            gGE0 = g >= 0;
            gLT0 = g < 0;
            uInf = ub == inf;
            lInf = lb == -inf;
            
            %calculate v including handling uncontrained elements
            v(gGE0) = x(gGE0) - lb(gGE0);
            v(gLT0) = x(gLT0) - ub(gLT0);
            v(gLT0 & uInf) = -1;
            v(gGE0 & lInf) = 1;
            
            %if v == 0, D becomes singular
            vAbs = abs(v);
            chk = vAbs < 100*eps; %not really sure on this choice
            vAbs(chk) = 100*eps;
            %D = diag(vAbs.^-0.5);
            D2 = diag(vAbs.^-1); %D*D
            
        end
        
        function [xHat,dfval,fvalHat] = nextStep(x,s,gamma,lb,ub,fval,H,f)
            %note dvfal < 0 means it has decreased (a good thing...)
            
            %create new step and ensure it is within bounds:
            xHat = x + gamma*s;
            xHat = qpBox.moveInBounds(xHat,lb,ub);
            %evaluate function value at new step and compare with previous
            fvalHat = qpBox.fQuad(H,f,xHat);
            dfval = fvalHat - fval;
        end
        
        function x = moveInBounds(x,lb,ub)
            %rounding errors might take xhat slightly outside box (can propagate and
            %blow up solution). Is there a way of scaling step so doesn't happen in first place?
            xGTub = x > ub;
            xLTlb = x < lb;
            x(xGTub) = ub(xGTub);
            x(xLTlb) = lb(xLTlb);
        end
    end
end