function [G, B] = toeplksolvetoeplk(G1, B1, G2, B2)
% function [G, B] = toeplksolvetoeplk(G1, B1, G2, B2)
%
% Solve a Toeplitz-like system with Toeplitz-like RHS, i.e. T1 \ T2.
%
% Input:
%   G1,B1 -- generator for T1
%   G2,B2 -- generator for T2
%
% Output:
%   G,B   -- generator for T1 \ T2
%
% Complexity is roughly O(d1*d2*n^2), where d1/d2 are the displacement
% ranks of the input.  This bound may be inaccurate, as it doesn't account
% for Gu's pivoting cost inside the GKO solve.

n = size(G1, 1);
e1 = zeros(n,1);
e1(1) = 1;

% Formulas based on Theorem X.X (FIXME) for generators of r(T).
G = [ratgen_solve(G1, B1, [-G1, G2]), e1];
Bpart = ratgen_solve(B1, G1, [B1, e1], B2, G2);
B = [Bpart(:,1:(end-1)), B2, Bpart(:,end)];

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


