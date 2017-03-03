classdef TLMatStein
    properties
        G
        B
    end
    
    methods
        function TL = TLMatStein(G,B)
            % TL = TLMatStein(A)   % Approximate A with TL matrix
            % TL = TLMatStein(c,r) % Generate TL form Toeplitz matrix
            % TL = TLMatStein(G,B) % Generators are directly given
        end
    end
end