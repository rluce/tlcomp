function E = tleye(n, nbis)
% E = tleye(n)  --  Identity matrix in TL format

if nargin == 1
    if ~isscalar(n)
        [mm,nn] = size(n);
        if mm~= 1
            error('tlcomp:InconsistentInput', 'Size vector must be row vector');
        end
        
        if nn > 2
            error('tlcomp:InconsistentInput', 'Size vector must have two elements');
        end
        
        if n(1)~=n(2)
            error('tlcomp:Unsupported', 'Only square matrices implemented');
        end
        
        n = n(1);
    end
elseif nargin == 2
    if ~isscalar(n) || ~isscalar(nbis)
        error('tlcomp:InconsistentInput', 'dimensions must be scalars');
    end
    
    if n ~= nbis
        error('tlcomp:Unsupported', 'Only square matrices implemented');
    end
else
    assert(false)
end
e1 = zeros(n,1);

if n > 0
    e1(1) = 1;
end

E = TLMat(e1,e1);
end