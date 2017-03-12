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
        
        function T = transpose(T)
            tmp = T.c;
            T.c = T.r;
            T.r = tmp;
        end
        
        function T = ctranspose(T)
            tmp = T.c;
            T.c = conj(T.r);
            T.r = conj(tmp);
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
                    elseif is_exact_toeplitz(op2)
                        % It's a dense toeplitz matrix
                        S = op1.add_dense_toeplitz_matrix(op2);
                    else
                        % It's some dense matrix
                        S = op1.add_dense_matrix(op2);
                    end
                case 'ToepMat'
                    S = op1.add_toepmat(op2);
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
        
        function S = add_toepmat(TM1, TM2)
            if any(size(TM1) ~= size(TM2))
                error('tlzstein:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            TM1.r = TM1.r + TM2.r;
            TM1.c = TM1.c + TM2.c;
            S = TM1;
        end
        
        function S = add_dense_toeplitz_matrix(TM, T)
            if any(size(TM) ~= size(T))
                error('tlzstein:InconsistentInput', ...
                'Matrix dimensions must agree');
            end

            if isempty(T)
                S = TM;
                return;
            end
                
            cc = T(:,1);
            rr = T(1,:);

            TM.c = TM.c + cc(:);
            TM.r = TM.r + rr(:);
            S = TM;
        end
        
        function S = add_dense_matrix(TM, A)
            if any(size(TM) ~= size(A))
                error('tlzstein:InconsistentInput', ...
                'Matrix dimensions must agree');
            end
            S = full(TM) + A;
        end
        
        
        function P = mtimes(op1, op2)
            if isa(op1, 'ToepMat')
                % Dispatch on first operand
                switch class(op2)
                    case 'double'
                        P = dispatch_tm_mtimes_double(op1, op2);
                    case 'TLMat'
                        % Promote op1, resort to TL matmul
                        TL = TLMat(op1.c, op1.r);
                        P = TL * op2;
                    case 'ToepMat'
                        % Promote op1, resort to TL matmul
                        TL = TLMat(op1.c, op1.r);
                        P = TL * op2;                        
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
            
            [m, n] = size(TM);
            if n~= mm
                error('tlzstein:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            
            if isempty(A)
                P = zeros(m, nn);
                return;
            end
            
            % A is at least 2x2
            P = toepmult(TM.c, TM.r, A);

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
            
            [m, n] = size(TM);
            if m~= nn
                error('tlzstein:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            
            if isempty(A)
                P = zeros(mm, n);
                return;
            end
            
            % A is at least 2x2
            TM = TM';
            P = toepmult(TM.c, TM.r, A')';

        end

        function P = mtimes_scalar(TM, s)
            TM.c = s * TM.c;
            TM.r = s * TM.r;
            P = TM;
        end
    end
end
