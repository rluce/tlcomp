function E = tleye(n)
% E = tleye(n)  --  Identity matrix in TL format

e1 = zeros(n,1);

if n > 0
    e1(1) = 1;
end

E = TLMat(e1,e1);
end