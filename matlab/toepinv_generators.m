function [Ginv, Binv] = toepinv_generators(c, r)
% [Ginv, Binv] = toepinv_generators(c, r)
%
% Computes a generator for the inverse of a Toeplitz matrix, with respect
% to the Sylvester displacement equation
%
%   Zp * T - T * Zm,
%
% with Zp/m being the unit +/-1 circulants.
%
% Input:
%   c,r -- first column/row of T
%
% Output:
%   Ginv, Binv -- generator for T^{-1}
%
% Cost: Two solves with each of T and T', currently done via drsolve,
% resulting in O(n^2).

% We use the following approach:  Consider the augmentation
%
%   M = [
%         -T, I;
%          I, 0
%       ]
%
% and the canoncial generator for T, G = [g, e1], B = [en, b].  It can be
% checked that GG = [g, e1; 2*e1, 0] and BB = [en, b; 0, 2*en] is a
% generator for M (w.r.t. Zp \oplus Zp / Zm \oplus Zm).  From the general
% Schur complement formula the generator expressions below then follow.

n = length(c);
assert(length(r) == n && c(1) == r(1));

% Match notation in the paper
[G, B] = toepgen(c,r);
e1 = G(:,2);
en = B(:,1);

Ginv = [ 2*e1, zeros(n,1) ] + toepsolve(c,r, -G);
Binv = [ zeros(n,1), -2*en ] + toepsolve(conj(r), conj(c), B);

end
