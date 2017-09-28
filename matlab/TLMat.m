classdef TLMat
    properties
        G
        B
    end
    
    methods
        function TL = TLMat(c, r, flag)
            % TL = TLMat(c ,r)
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
                            error('tlzstein:InconsistentInput', ...
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
                        error('tlzstein:InconsistentInput', ...
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
            tmp = TL.G;
            TL.G = conj(TL.B);
            TL.B = conj(tmp);
        end

        function TL = ctranspose(TL)
            tmp = TL.G;
            TL.G = TL.B;
            TL.B = tmp;
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
                    error('tlzstein:NotImplemented', ...
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
                error('tlzstein:InconsistentInput', ...
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
                error('tlzstein:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            
            B = toeplkmult(TL.G, TL.B, A);
        end
        
        function P = mtimes_tlmat(TL1, TL2)
            if any(size(TL1) ~= size(TL2))
                error('tlzstein:InconsistentInput', ...
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
                error('tlzstein:InconsistentInput', ...
                    'Matrix dimensions must agree');
            end
            
            switch class(op)
                case 'double'
                    D = toeplksolve(TL.G, TL.B, op);
                case 'TLMat'
                    D = TL.mldivide_tlmat(op);
                case 'ToepMat'
                    op = TLMat(op.c, op.r);
                    D = TL.mldivide(op);
                otherwise
                    error('tlzstein:NotImplemented', ...
                        'No mldivide for this operand type, fixme');
            end
        end
        
        function D = mldivide_tlmat(TL, TL_rhs)
            [Gs, Bs] = toeplksolvetoeplk(TL.G, TL.B, TL_rhs.G, TL_rhs.B);
            D = TLMat(Gs, Bs, 'GB');
        end
        
        
        function TL = compress(TL)
            [TL.G, TL.B] = gencompress(TL.G, TL.B);
        end

        %%%%%%%% CAUTION DEBUG ONLY
        function d = det(TL)
            warning('tlzstein:CubicOperatoin', ...
                'det not implemented, using dense det');
            d = det(full(TL));
        end
        
        function TL = mrdivide(TL, s)
           % Implement only a special case for debug purpose
           assert(isa(TL, 'TLMat'));
           assert(isa(s, 'double'));
           assert(all(size(s) == [1,1]));
           TL = scalar_mult(TL, 1.0/s);
        end
        
        function val = norm(TL, p)
            val = norm(full(TL), p);
        end
    
        function TL = inv(TL)
            [m,n] = size(TL);
            TL = TL \ tleye(n);
        end
    end % of methods section
end
