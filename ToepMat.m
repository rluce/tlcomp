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
                case 'double'
                    [m,n] = size(op2);
                    if m==1 && n==1
                        % Add a scalar to all entries
                        S = op1.add_scalar(op2);
                    else
                        % It's a matrix
                        S = op1.add_dense_matrix(op2);
                    end
                otherwise
                    error('tlzstein:NotImplemented', ...
                        'Addition not implemented for this operand');
            end
        end
        
        function S = add_scalar(TM, alpha)
            TM.r = TM.r + alpha;
            TM.c = TM.c + alpha;
            S = TM;
        end
        
        function P = mtimes(op1, op2)
            if isa(op1, 'ToepMat')
                % Dispatch on first operand
                switch class(op2)
                    case 'double'
                        P = dispatch_tm_mtimes_double(op1, op2);
                    otherwise
                        error('tlzstein:NotImplemented', ...
                            'Multiplication not implemented for this operand');
                end
            elseif isa(op2, 'ToepMat')
                % Dispatch on second operand
                switch class(op1)
                    case 'double'
                        P = dispatch_double_mtimes_tm(op1, op2);
                    otherwise
                        error('tlzstein:NotImplemented', ...
                            'Multiplication not implemented for this operand');

                end
            else
                % Impossible
                assert(false);
            end
        end
        
        function P = dispatch_tm_mtimes_double(TM, A)
            assert(isa(TM, 'ToepMat'));
            assert(isa(A, 'double'));
            
            [mm,nn] = size(A);
            if mm == 1 && nn == 1
                % Scalar multiplication, which is commutative
                P = TM.mtimes_scalar(A);
                return;
            end
        end

        function P = dispatch_double_mtimes_tm(A, TM)
            assert(isa(TM, 'ToepMat'));
            assert(isa(A, 'double'));
            
            [mm,nn] = size(A);
            if mm == 1 && nn == 1
                % Scalar multiplication, which is commutative
                P = TM.mtimes_scalar(A);
                return;
            end
        end

        
        function P = mtimes_scalar(TM, s)
            TM.c = s * TM.c;
            TM.r = s * TM.r;
            P = TM;
        end
    end
end
