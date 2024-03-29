classdef ToepMat
    properties
        c % First column
        r % First row
    end
    
    methods
        function T = ToepMat(c, r)
            % T = ToepMat(c,r)
            %   T is Toeplitz matrix with prescribed first col/row
            if nargin == 2
                if min(size(c)) > 1 || min(size(r)) > 1
                    error('tlcomp:InconsistentInput', 'c, r must be vectors');
                end
                
                % Homogenize
                c = c(:);
                r = r(:);
                
                if length(c) ~= length(r)
                    error('tlcomp:InconsistentInput', 'c, r must have same length');
                end
                
                n = length(r);
                
                if n>0 && c(1) ~= r(1)
                    warning('tlcomp:InconsistentInput', ...
                        'First elements of c,r are not equal, c wins conflict');
                    r(1) = c(1);
                end
                
                T.c = c;
                T.r = r;
            else
                % Only r is provided
                assert(nargin == 1);
                r = c;
                if min(size(r)) > 1
                    error('tlcomp:InconsistentInput', 'r must be a vector');
                end
                
                if isreal(r(1))
                    c = conj(r);
                else
                    c = [r(1); conj(r(2:end))];
                end
                T.c = c(:);
                T.r = r(:);
            end
        end
        
        function [s1, s2] = size(T, dim)
            if nargin == 1
                s1 = length(T.c);
                s2 = length(T.r);
                if nargout < 2
                    s1 = [s1, s2];
                end
            elseif nargin == 2
                if dim == 1
                    s1 = length(T.c);
                elseif dim == 2
                    s1 = length(T.r);
                else
                    error('tlcomp:InconsistentInput', ...
                        'No such dimension');
                end
            else
                % Cannot happen unless code broken
                assert(false);
            end
            
        end

        function disp(T)
            fprintf('  %dx%d ToepMat\n', size(T, 1), size(T, 2));
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
                    elseif istoeplitz(op2)
                        % It's a dense toeplitz matrix
                        S = op1.add_dense_toeplitz_matrix(op2);
                    else
                        % It's some dense matrix
                        S = op1.add_dense_matrix(op2);
                    end
                case 'ToepMat'
                    S = op1.add_toepmat(op2);
                otherwise
                    error('tlcomp:NotImplemented', ...
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
                error('tlcomp:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            TM1.r = TM1.r + TM2.r;
            TM1.c = TM1.c + TM2.c;
            S = TM1;
        end
        
        function S = add_dense_toeplitz_matrix(TM, T)
            if any(size(TM) ~= size(T))
                error('tlcomp:InconsistentInput', ...
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
                error('tlcomp:InconsistentInput', ...
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
                        error('tlcomp:NotImplemented', ...
                            'Multiplication not implemented for this operand');
                end
            elseif isa(op2, 'ToepMat')
                % Dispatch on second operand
                switch class(op1)
                    case 'double'
                        P = dispatch_double_mtimes_tm(op1, op2);
                    otherwise
                        error('tlcomp:NotImplemented', ...
                            'Multiplication not implemented for this operand');
                        
                end
            else
                % Impossible
                assert(false);
            end
        end
        
        
        function R = chol(TM)
            % Simplistic Schur algorithm for Cholesky factorization
            % DEBUG ONLY
            
            n = size(TM,1);
            
            % Normalization, not really needed but handy for data
            % inspection later.
            t0 = TM.c(1);
            G = 1./t0 * [TM.r, [0;TM.r(2:end)]]';
            
            R = zeros(n,n);
            
            % Init: Read off first row of chol(T), prepare Schur complement
            % generator
            R(1,1:n) = G(1,1:n);
            G(1,2:n) = G(1,1:n-1);
            G(1,1) = 0;
            
            % Loop to generate rows 2,...,n of R
            for i = 2:n
                % Parameter of the rotation
                rho = -G(2,i) / G(1,i);
                s = 1.0 / sqrt((1-rho)*(1+rho));
                
                % Rotation / elimination step
                G(:, i:n) = s * [1, rho; rho, 1] * G(:, i:n);
                
                % Copy over first row of current Schur complement
                R(i, i:n) = G(1,i:n);
                
                % Shift generator
                G(1,i+1:n) = G(1,i:n-1);
                G(1,i) = 0; % not really needed but nice for inspection
            end
            
            R = sqrt(t0) * R;
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
                error('tlcomp:InconsistentInput', ...
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
                error('tlcomp:InconsistentInput', ...
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
        
        function P = polyvalm(p, TM)
            [GP, BP] = toeppolyvalm(TM.c, TM.r, p);
            P = TLMat(GP, BP, 'GB');
        end
        
        function TMinv = inv(TM)
            [Ginv, Binv] = toepinvgen(TM.c, TM.r);
            TMinv = TLMat(Ginv, Binv, 'GB');
        end
        
        function R = ratevalm(TM, p, q)
            [Gr, Br] = toepratvalm(TM.c,TM.r,p,q);
            R = TLMat(Gr, Br, 'GB');
        end
        
        function X = mldivide(TM, B)
            assert(isa(TM, 'ToepMat'));
            
            [m,~] = size(TM);
            [mm, ~] = size(B);
            if m ~= mm
                error('tlcomp:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            
            switch class(B)
                case 'double'
                    X = toepsolve(TM.c, TM.r, B);
                case 'ToepMat'
                    % Promote to TLMat, handle solve there
                    TL = TLMat(TM.c, TM.r);
                    X = TL \ B;
                case 'TLMat'
                    % Promote to TLMat, handle solve there
                    TL = TLMat(TM.c, TM.r);
                    X = TL \ B;
                otherwise
                    error('tlcomp:NotImplemented', ...
                        'Operation not implemented yet, fixme');
            end
        end
        
        function val = norm(TM, p)
            % v = norm(T, p) -- computes the matrix p-norm
            %
            % Supported values for p are 1, inf and 'fro'.
            %
            % Warning: if no value for p is given, or p equals to 2, the
            % entailing computation is carried out with Matlab's unstructured
            % matrices, which takes a cubic number of operations.  Consider
            % using normest instead.
            %
            % See also: normest.m norm.m
            
            if nargin < 2
                p = 2;
            end
            
            if p == 2
                warning('tlcomp:Unsupported', ...
                    ['2-norm not supported, computation will be slow. ', ...
                    'Consider using normest']);
            end

            switch p
                case 1
                    val = toepnorm(TM.c, TM.r, p);
                case 2
                    val = norm(full(TM), 2);
                case inf
                    val = toepnorm(TM.c, TM.r, p);
                case 'inf'
                    val = toepnorm(TM.c, TM.r, inf);
                case 'fro'
                    val = toepnorm(TM.c, TM.r, 'fro');
                otherwise
                    error('tlcomp:InconsistentInput', 'Unsupported matrix norm');
            end
        end
        
        function l = length(TM)
            l = max(size(TM));
        end
        
        function R = mrdivide(op1, op2)
            if ~isa(op1, 'ToepMat')
                error('tlcomp:NotImplemented', 'Not implemented, fixme.');
            end
            
            switch class(op2)
                case 'double'
                    R = dispatch_mrdivide_tm_double(op1, op2);
                case 'TLMat'
                    % Promote self to TLMat, handle mldivide there
                    TL1 = TLMat(op1.c, op1.r);
                    R = TL1 / op2;
                otherwise
                    error('tlcomp:NotImplemented', 'Not implemented, fixme.');
            end
        end
        
        function R = dispatch_mrdivide_tm_double(TM, A)
            assert(isa(TM, 'ToepMat'));
            assert(isa(A, 'double'));
            
            [mm, nn] = size(A);
            if mm == 1 && nn == 1
                % Resolves through scalar multiplication
                R = TM.mtimes_scalar(1.0 / A);
                return;
            end
            error('tlcomp:NotImplemented', 'Not implemented, fixme.');
        end
    
        function Ts = mpower(base, exponent)
            if ~isa(base, 'ToepMat')
                error('tlcomp:NotImplemented', 'Not implemented, fixme.');
            end
            
            switch class(exponent)
                case 'double'
                    if (exponent - floor(exponent) == 0) && (exponent >= 0)
                        Ts = mpower_TM_integer(base, exponent);
                    else
                        error('tlcomp:NotImplemented', 'Not implemented, fixme.');
                    end                
                otherwise
                    error('tlcomp:NotImplemented', 'Not implemented, fixme.');
            end
        end
        
        function Ts = mpower_TM_integer(TM, s)
            if s==0
                Ts = toepeye(size(TM));
            elseif s==1
                Ts = TM;
            else
                [Gpow, Bpow] = toeppowers(TM.c, TM.r, s);
                Ts = TLMat(Gpow{end}, Bpow{end}, 'GB');
            end
        end
        
        function val = det(TM)
            [G, B] = toepgen(TM.c, TM.r);
            TL = TLMat(G, B, 'GB');
            val = det(TL);
        end
    end
end
