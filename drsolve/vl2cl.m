function [G,H,x,b] = vl2cl(x,eta,G,H)
%VL2CL	Converts a square Vandermonde-like matrix to Cauchy-like.
%   [G,H,a,b] = V2CL(x,eta,Gv,Hv) returns the quantities used in the
%   displacement equation
%
%       diag(a) * C - C * diag(b) = G * H'
%
%   of a Cauchy-like matrix C corresponding to the Vandermonde-like
%   matrix V with displacement equation
%
%       diag(x) * V - V * Z(eta)' = Gv * Hv'.
%
%   Here Z(eta) is the square eta-cyclic forward shift matrix defined as
%   Z(eta) = [O,eta;eye(n-1),0], where 0 are null blocks.
%
%   The matrices C and V are such that C = V*F(eta), where
%   F(eta)=FTIMES(n,[],eta) is a unitary matrix.
%
%   See also v2cl, ftimes.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

if (nargin~=4), error('drsolve:v2cl:nargin','wrong arguments number'); end
x = x(:);
n = size(x,1);
b = conj(nroots1(n,eta));
if size(intersect(b,x))>1
    error('drsolve:vl2cl','C non ricostruibile, cambia eta');
end
H = ftimes(H,'A',eta);
