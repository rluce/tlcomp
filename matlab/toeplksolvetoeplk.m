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
% Complexity is O(d1*d2*n^2), where d1/d2 are the displacement
% ranks of the input matrices

n = size(G1, 1);
e1 = zeros(n,1);
e1(1) = 1;
en = zeros(n,1);
en(n) = 1;

% Formulas based on Theorem X.X (FIXME) for generators of r(T).
G = [ toeplksolve(G1, B1, [G1, G2]), -2 * e1 ];

Bpart = toeplksolve(G1, B1, [B1, en], true);
Bpart = toeplkmult(G2, B2, Bpart, true);
B = [-Bpart(:,1:(end-1)), B2, -Bpart(:,end)];

end

