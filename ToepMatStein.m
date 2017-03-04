classdef ToepMatStein
    properties
        c % First column
        r % First row
    end
    
    methods
        function T = ToepMatStein(c, r)
            % T = ToepMatStein(c,r)
            %   T is Toeplitz matrix with prescribed first col/row
            
            
            if min(size(c)) > 1 || min(size(r)) > 1
                error('funmd:InconsistentInput', 'c, r must be vectors');
            end
            
            % Homogenize
            c = c(:);
            r = r(:);
            
            if length(c) ~= length(r)
                error('funmd:InconsistentInput', 'c, r must have same length');
            end
            
            n = length(r);
            
            if n>0 && c(1) ~= r(1)
                error('funmd:InconsistentInput', 'c(1),r(1) must be equal');
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
        
    end
end