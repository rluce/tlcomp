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
                [G, B] = compute_generator(c);
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
                        % We assume it to be first row/col
                        [G, B] = stein_generator(c,r);
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
        
        function [m,n] = size(TL)
            m = size(TL.G, 1);
            n = size(TL.B, 1);
            if nargout == 1
                m = [m,n];
            end
        end
        
        function d = drank(TL)
            d = size(TL.G, 2);
        end
        
        function T = full(TL)
            T = stein_reconstruction(TL.G, TL.B);
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
            
            % Structured + unstructured = unstructured, in the absence
            % of further information.  Resort to full matrices.
            S = full(TL) + A;
        end
        
        function S = add_tlmat(T1, T2)
            assert(isa(T1, 'TLMat') && isa(T2, 'TLMat'));

            T1.G = [T1.G, T2.G];
            T1.B = [T1.B, T2.B];

            S = T1;
            S = S.compress;
        end

        function S = uminus(TL)
            TL.B = -TL.B;
            S = TL;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function TL = compress(TL)
            [TL.G, TL.B] = gencompress(TL.G, TL.B);
        end
        
    end % of methods section
end
