classdef ToepMat
    properties
        c % First column
        r % First row
    end
    
    methods
        function T = ToepMat(c, r)
            % T = ToepMat(c,r)
            %   T is Toeplitz matrix with prescribed first col/row
            
            if min(size(c)) > 1 || min(size(r)) > 1
                error('tlzstein:InconsistentInput', 'c, r must be vectors');
            end
            
            % Homogenize
            c = c(:);
            r = r(:);
            
            if length(c) ~= length(r)
                error('tlzstein:InconsistentInput', 'c, r must have same length');
            end
            
            n = length(r);
            
            if n>0 && c(1) ~= r(1)
                error('tlzstein:InconsistentInput', 'c(1),r(1) must be equal');
            end
            
            T.c = c;
            T.r = r;
        end
       
        function [s1, s2] = size(T)
            s1 = length(T.c);
            s2 = length(T.r);
            if nargout < 2
                s1 = [s1, s2];
            end
        end
        
        function TF = full(T)
            if isempty(T.c)
                TF = [];
                return;
            end
            
            TF = toeplitz(T.c, T.r);
        end
        
        function T = uplus(T)
            % Do nothing
        end
        
        function T = uminus(T)
            T.c = -T.c;
            T.r = -T.r;
        end
        
        function S = minus(op1, op2)
            S = op1 + (-op2);
        end
        
        function S = plus(op1, op2)
            if ~isa(op1, 'ToepMat')
                % Plus is commutative, make ToepMat first operand
                tmp = op1;
                op1 = op2;
                op2 = tmp;
            end
            
            assert(isa(op1, 'ToepMat'));
            % Dispatch
            switch class(op2)
                case 'TLMat'
                    % Promote op1 to TLMat, and add
                    S = TLMat(op1.c, op1.r) + op2;
                otherwise
                    error('tlzstein:NotImplemented', ...
                        'Addition not implemented for this operand');
            end
        end
        
    end
end
