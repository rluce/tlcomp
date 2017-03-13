function [G, B] = toepratvalm(c, r, p, q)
% [G, B] = toepratvalm(c,r,p,q)
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
% and that r(z) has no poles at or near the eigengalues of T.

n = length(c);
e1 = zeros(n,1);
e1(1) = 1;

% TODO
% This is inefficient, as one single et of polynomials needs only to be
% computed.
[Gp, Bp] = toeppolyvalm(c, r, p);
[Gq, Bq] = toeppolyvalm(c, r, q);

% Compress as much as we can
[Gp, Bp] = gencompress(Gp, Bp);
[Gq, Bq] = gencompress(Gq, Bq);

% Formulas based on Theorem X.X (FIXME) for generators of r(T).
G = [ratgen_solve(Gq, Bq, [-Gq, Gp]), e1];
Bpart = ratgen_solve(Bq, Gq, [Bq, e1], Bp, Gp);
B = [Bpart(:,1:(end-1)), Bp, Bpart(:,end)];

end


% Utiility function for rat generator construction
function x = ratgen_solve(G, B, rhs, multG, multB)
% solve TL system with G,B generated matrix and given rhs.  Optionally, do
% an intermediate TL multiplication.

x = vapply(rhs, 'inv');
x = toeplksolve(G, B, x);
if nargin >= 5
    x = toeplkmult(multG, multB, x);
end

x = vapply(x);

end


