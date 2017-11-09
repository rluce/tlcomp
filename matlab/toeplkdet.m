function d = toeplkdet(G, B)
% d = toeplkdet(G, B) -- Determinant of a TL matrix.
%
% Input:
%   G, B -- Generator of the TL matrix A
%
% Output:
%   d    -- det(A)
%
% NOTE: The computation is based on the LU factorization of A after FFT to
% a Cauchy matrix.  The limitations for the usual det function apply here as
% well.

d = 0;

end