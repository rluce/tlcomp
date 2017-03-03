classdef TLMatStein
    properties
        G
        B
    end
    
    methods
        function TL = TLMatStein(c, r, flag)
            % TL = TLMatStein(c ,r)
            %   Generate a TL matrix from the first col/row of a Toeplitz
            %   matrix.
            %
            % TL = TLMatStein(G, B)
            %   Generate a TL matrix with prescribed generators of size
            %   n-by-r.
            %
            % TL = TLMatStein(G, B, 'GB')
            %   Force intepretation as generators (only useful if r=1).
            %
            % TL = TLMatStein(A)
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
                            error('funmd:InconsistentInput', ...
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
                        error('funmd:InconsistentInput', ...
                            'Generator matrices must be of equal size');
                    end
                    G = c;
                    B = r;
                end
            end
            TL.G = G;
            TL.B = B;
        end % of constructor
        
        function [m,n] = size(TL)
            m = size(TL.G, 1);
            n = size(TL.B, 1);
            if nargout == 1
                m = [m,n];
            end
        end
    end % of methods section
end