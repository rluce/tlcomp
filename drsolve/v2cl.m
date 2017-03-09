function [G,H,x,b] = v2cl(x,eta)
%V2CL	Converts a square Vandermonde matrix to Cauchy-like.
%   [G,H,a,b] = V2CL(x,eta) returns the quantities used in the
%   displacement equation
%
%       diag(a) * C - C * diag(b) = G * H'
%
%   of a Cauchy-like matrix C corresponding to the Vandermonde matrix
%   V=VANDER(x).
%   The matrices C and V are such that C = V*F(eta), where
%   F(eta)=FTIMES(n,[],eta) is a unitary matrix.
%
%   See also vl2cl, ftimes, vander.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

if (nargin<2) || isempty(eta),  eta = 1;  end
x = x(:);
n = size(x,1);

b = conj(nroots1(n,eta));
if numel(intersect(b,x))>0
    error('drsolve:v2cl','C Non ricostruibile cambia eta');
end

G = x.^n-conj(eta);
H = repmat(n^(-1/2),[n 1]);
