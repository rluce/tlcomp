function [G,H,a,b] = thl2cl(G,H)
%THL2CL	Converts a square Toeplitz+Hankel-like matrix to Cauchy-like.
%   [G,H,a,b] = T2CL(G_THL,H_THL) returns the quantities used in the
%   displacement equation
%
%       diag(a) * C - C * diag(b) = G * H'
%
%   of a Cauchy-like matrix C corresponding to the
%   Toeplitz+Hankel-like matrix M. The matrices C and M are such that
%   C = U*M*V, where U=STIMES(n) and V=CTIMES(n) are unitary matrices.
%
%   See also th2cl, ctimes, stimes.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

if nargin~=2, error('drsolve:thl2cl:nargin','wrong arguments number'), end
if any(size(G)~=size(H)), error('drsolve:thl2cl:nargin','wrong argument size'); end

n = size(G,1);

a = 2*cos((1:n).'*(pi/(n+1)));
b = 2*cos((0:n-1).'*(pi/n));
G = stimes(G);% ~dst
H = ctimes(H,'T');% dct
