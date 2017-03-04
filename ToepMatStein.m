classdef ToepMatStein
    properties
        c % First column
        r % First row
    end
    
    methods
        function T = ToepMatStein(c, r)
            % T = ToepMatStein(c,r)
            %   T is Toeplitz matrix with prescribed first col/row
            
            T.c = c;
            T.r = r;
        end
        
    end
end