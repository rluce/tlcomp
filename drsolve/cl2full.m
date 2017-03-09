function [C] = cl2full(G,H,t,s)
%CL2FULL Construct a Cauchy-like matrix from its generators.
%   C = CL2FULL(G,H,t,s) computes the Cauchy-like matrix C from its
%   displacement equation
%
%      diag(t) * C - C * diag(s) = G * H'.
%
%   The matrix C must be fully-reconstructable, i.e., t(i)<>s(j) for
%   any (i,j).
%
%   See also cltimes.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

t = t(:);
s = s(:);

C = intersect(t,s);
if numel(C)
    warning('drsolve:cl2full', 'C is not fully reconstructable');
end

[T,S] = ndgrid(t,s);
C = (G*H') ./ (T-S);

