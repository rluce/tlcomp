function [G, B] = toepratvalm(c, r, p, q)
% [G, B] = toepratvalm(c, r, p, q)
%
% Evaluate rational function of a Toeplitz matrix, i.e. r(T) =p(T) / q(T).
%
% Input:
%   c,r   -- first column and row of the Toeplitz matrix T
%   p,q   -- numerator and denominator polynomial
%
% Output:
%   G,B   -- generator for p(T)/q(T)
%
% NOTE:  It is assumed that r(z) = p(z)/q(z) is proper (no common roots)
% and that r(z) has no poles at or near the eigenvalues of T.

% Evaluate polynomials individually
[Gp, Bp] = toeppolyvalm(c, r, p);
[Gq, Bq] = toeppolyvalm(c, r, q);

% Compress as much as we can
[Gp, Bp] = gencompress(Gp, Bp);
[Gq, Bq] = gencompress(Gq, Bq);

[G, B] = toeplksolvetoeplk(Gq, Bq, Gp, Bp);

end
