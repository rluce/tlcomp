classdef TLMat
    properties
        G
        B
    end
    
    methods
        function TL = TLMat(c, r, flag)
            % TL = TLMat(c, r)
            %   Generate a TL matrix from the first col/row of a Toeplitz
            %   matrix.
            %
            % TL = TLMat(G, B)
            %   Generate a TL matrix with prescribed generators of size
            %   n-by-r.
            %
            % TL = TLMat(G, B, 'GB')
            %   Force intepretation as generators (only useful if r=1).
            %
            % TL = TLMat(A)
            %   Approximate a given matrix A with a Toeplitz matrix.
            %   (EXPENSIVE, for debug only).
            
            if nargin > 2 && strcmp(flag, 'GB')
                force_gb = true;
            else
                force_gb = false;
            end
            
            if nargin == 1
                [G, B] = toeplkgen(c);
            else
                % Two argument instantiation
                sc = size(c);
                sr = size(r);
                
                if min(sc) == 1 && min(sr) == 1
                    % Input was two vectors
                    if force_gb
                        % But we read it as generators
                        if sc(2) ~= sr(2)
                            error('tlcomp:InconsistentInput', ...
                                'Generator matrices must have equal length');
                        end
                        G = c;
                        B = r;
                    else
                        % We assume it to be first col/row
                        [G, B] = toepgen(c, r);
                    end
                else
                    % Two inputs, and at least one is not a vector,
                    % so the input is a generator.
                    if sc(2) ~= sr(2)
                        error('tlcomp:InconsistentInput', ...
                            'Generator matrices must be of equal size');
                    end
                    G = c;
                    B = r;
                end
            end
            TL.G = G;
            TL.B = B;
            TL = TL.compress();
        end % of constructor
        
        function [m,n] = size(TL, dim)
            if nargin == 1
                m = size(TL.G, 1);
                n = size(TL.B, 1);
            
                if nargout == 1
                    m = [m,n];
                end
                return;
            elseif nargin == 2
                if dim == 1
                    m = size(TL.G, 1);
                elseif dim == 2
                    m = size(TL.B, 1);
                else
                    assert(false); % Should raise an error....
                end
                return;
            end
            
        end
        
        function d = drank(TL)
            d = size(TL.G, 2);
        end
        
        function T = full(TL)
            T = toeplkreconstruct(TL.G, TL.B);
        end
        
        function TL = transpose(TL)
            TL = TL';
            TL.G = conj(TL.G);
            TL.B = conj(TL.B);
        end

        function TL = ctranspose(TL)
            [Gt, Bt] = toeplkctranspose(TL.G, TL.B);
            TL.G = Gt;
            TL.B = Bt;
            % Generator rank may have increased, but rank is in fact the
            % same.
            TL = TL.compress();
        end

        
        function S = minus(op1, op2)
            S = op1 + (-op2);
        end
        
        function S = plus(op1, op2)
            % Homogenize
            if ~isa(op1, 'TLMat')
                tmp = op2;
                op2 = op1;
                op1 = tmp;
            end
            
            % Dispatch
            switch class(op2)
                case 'TLMat'
                    S = op1.add_tlmat(op2);
                case 'double'
                    S = op1.add_double(op2);
                case 'ToepMat'
                    S = op1.add_toepmat(op2);
                otherwise
                    error('tlcomp:NotImplemented', ...
                        'Addition not implemented for this operand');
            end
        end
        
        function S = add_double(TL, A)
            assert(isa(A, 'double'));
            [m,n] = size(A);
            if m==1 && n==1
                if A==0
                    S = TL;
                    return;
                end
                % Interpret as entrywise addition
                augvec = A * ones(size(TL.B, 1), 1);
                A = TLMat(augvec, augvec);
                S = TL.add_tlmat(A);
                return;
            end
                
            % A is not scalar.
            [mm, nn] = size(TL);
            if mm~=m || nn~=n
                error('tlcomp:InconsistentInput', ...
                    'Matrix dimensions must agree.');
            end
            
            if n==0 && m ==0
                % Empty matrix
                S = TL;
                return;
            end
            
            % At this point A is a non-scalar, non-empty matrix
            
            if istoeplitz(A)
                % Magic: If A is a dense toeplitz matrix, we retain TL
                % structure
                S = TL.add_tlmat(TLMat(A(:,1), A(1,:)));
            else
                % Structured + unstructured = unstructured, in the absence
                % of further information.  Resort to full matrices.
                S = full(TL) + A;
            end
        end
        
        function S = add_tlmat(T1, T2)
            assert(isa(T1, 'TLMat') && isa(T2, 'TLMat'));

            T1.G = [T1.G, T2.G];
            T1.B = [T1.B, T2.B];

            S = T1;
            S = S.compress;
        end
        
        function S = add_toepmat(TL, TM)
            S = TL.add_tlmat(TLMat(TM.c, TM.r));
        end

        function S = uminus(TL)
            TL.B = -TL.B;
            S = TL;
        end
        
        function P = mtimes(op1, op2)
            if isa(op1, 'TLMat')
                switch class(op2)
                    case 'double'
                        [m,n] = size(op2);
                        if m==1 && n == 1
                            P = scalar_mult(op1, op2);
                        else
                            P = mtimes_double(op1, op2);
                        end
                    case 'TLMat'
                        P = op1.mtimes_tlmat(op2);
                    case 'ToepMat'
                        TL = TLMat(op2.c, op2.r);
                        P = op1.mtimes_tlmat(TL);
                    otherwise
                        error('Not implemented');
                end
            elseif isa(op2, 'TLMat')
                switch class(op1)
                    case 'double'
                        [m,n] = size(op1);
                        if m==1 && n == 1
                            P = scalar_mult(op2, op1);
                        else
                            P = mtimes_double(op2', op1')';
                        end

                    otherwise
                        error('Not implemented');
                end
            else
                assert(false);
            end
        end
        
        function TL = scalar_mult(TL, s)
            TL.B = conj(s) * TL.B;
        end
        
        function B = mtimes_double(TL, A)
            if size(TL,2) ~= size(A,1)
                error('tlcomp:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            
            B = toeplkmult(TL.G, TL.B, A);
        end
        
        function P = mtimes_tlmat(TL1, TL2)
            if any(size(TL1) ~= size(TL2))
                error('tlcomp:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            [Gp, Bp] = toeplkprod(TL1.G, TL1.B, TL2.G, TL2.B);
            P = TLMat(Gp, Bp, 'GB');
        end
        
        function D = mldivide(TL, op)
            % Should only be called with first operand TL?
            % If this ever triggers, the dispatch must be extended.
            assert(isa(TL, 'TLMat'));
           
            m = size(TL, 1);
            mm = size(op,1);
            if mm ~= m
                error('tlcomp:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            
            switch class(op)
                case 'double'
                    D = toeplksolve(TL.G, TL.B, op);
                case 'TLMat'
                    D = TL.mldivide_tlmat_tlmat(op);
                case 'ToepMat'
                    op = TLMat(op.c, op.r);
                    D = TL.mldivide(op);
                otherwise
                    error('tlcomp:NotImplemented', ...
                        'No mldivide for this operand type, fixme!');
            end
        end
        
        function D = mldivide_tlmat_tlmat(TL, TL_rhs)
            [Gs, Bs] = toeplksolvetoeplk(TL.G, TL.B, TL_rhs.G, TL_rhs.B);
            D = TLMat(Gs, Bs, 'GB');
        end
        
        
        function TL = compress(TL)
            [TL.G, TL.B] = gencompress(TL.G, TL.B);
        end

        function result = mrdivide(op_left, op_right)
            % This function only dispatches and (maybe) verifies input.  No
            % computation here.
            
            if isa(op_left, 'TLMat')
                switch class(op_right)
                    case 'double'
                        if isscalar(op_right)
                            result = mrdivide_TL_scalar(op_left, op_right);
                        else
                            assert(false, 'Code path not exising, fixme');
                        end
                    case 'TLMat'
                        result = mrdivide_TL_TL(op_left, op_right);
                    otherwise
                        error('tlcomp:NotImplemented', ...
                            'No mrdivide for this operand type, fixme!');
                        
                end
            else
                assert(isa(op_right, 'TLMat'));
                
                switch class(op_left)
                    case 'double'
                        result = mrdivide_double_TL(op_left, op_right);
                    otherwise
                        error('tlcomp:NotImplemented', ...
                            'No mrdivide for this operand type, fixme!');
                end

            end
            
        end
        
        function TL = mrdivide_TL_scalar(TL, s)
            % Compute TL / (double scalar)
            TL = scalar_mult(TL, 1.0/s);
        end

        function result = mrdivide_double_TL(A, TL)
            % Compute (double matrix) / TL, result is a double matrix
            if size(A, 2) ~= size(TL, 2)
                error('tlcomp:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            
            result = toeplksolve(TL.G, TL.B, A', true)';
        end
        
        function result = mrdivide_TL_TL(TL1, TL2)
            % Compute TL / TL, result is a TL
            % We do it the naive way TL1/TL2 = (TL2'\TL1')'
            %
            % TODO this could easily made better:  extend toeplksolvetoeplk
            % to accept ctransposes of either argument, and thus switch
            % implicitly between Zp/Zm and Zm/Zp displacement.  This avoids
            % the numerical(!) ctranspose ops that are needed here!
            %
            % This should definitly improved, see issue #7.
            result = mldivide(TL2', TL1')';
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
                    val = toeplknorm(TM.G, TM.B, p);
                case 2
                    val = norm(full(TM), 2);
                case inf
                    val = toeplknorm(TM.G, TM.B, p);
                case 'inf'
                    val = toeplknorm(TM.G, TM.B, inf);
                case 'fro'
                    val = toeplknorm(TM.G, TM.B, 'fro');
                otherwise
                    error('tlcomp:InconsistentInput', 'Unsupported matrix norm');
            end
        end

        function Ts = mpower(base, exponent)
            if ~isa(base, 'TLMat')
                error('tlcomp:NotImplemented', 'Not implemented, fixme.');
            end
            
            switch class(exponent)
                case 'double'
                    if (exponent - floor(exponent) == 0) && (exponent >= 0)
                        Ts = mpower_TL_integer(base, exponent);
                    else
                        error('tlcomp:NotImplemented', 'Not implemented, fixme.');
                    end                
                otherwise
                    error('tlcomp:NotImplemented', 'Not implemented, fixme.');
            end
        end
        
        function Ts = mpower_TL_integer(TL, s)
            if s==0
                Ts = tleye(size(TL));
            else
                Ts = TL;
                for k=2:s
                    Ts = Ts * TL;
                end
            end
        end
        
        function val = gennorm(TL)
            [~, RG] = qr(TL.G,0);
            [~, RB] = qr(TL.B,0);
            val = norm(RG*RB');
        end
        
        %%%%%%%% CAUTION DEBUG ONLY
        function d = det(TL)
            warning('tlcomp:CubicOperation', ...
                'det not implemented, using dense det');
            d = det(full(TL));
        end
        
        function TL = inv(TL)
            [~,n] = size(TL);
            TL = TL \ tleye(n);
        end
    end % of methods section
end
