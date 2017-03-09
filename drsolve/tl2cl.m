function [GC,HC,a,b] = tl2cl(GT,HT,xi,eta)
%TL2CL	Converts a Toeplitz-like matrix to Cauchy-like.
%   [GC,HC,a,b] = TL2CL(GT,HT,xi,eta) returns the quantities used in
%   the displacement equation
%
%	   diag(a) * C - C * diag(b) = GC * HC'
%
%   of a Cauchy-like matrix C corresponding to the Toeplitz-like
%   matrix T whose displacement equation is
%
%      Zn(xi) * T - T * Zn(eta) = GT * HT',
%
%   Here Z(xi) is the square xi-cyclic forward shift matrix defined as
%   Z(xi) = [O,xi;eye(m-1),0], where 0 are null blocks.
%   The matrices C and T are such that C = F(xi)'*T*F(eta), where
%   F(xi)=FTIMES(m,[],xi) and F(eta)=FTIMES(n,[],eta) are unitary
%   matrices, and [m,n]=size(T).
%
%   See also t2tl, t2cl, ftimes.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

m  = size(GT,1);
n  = size(HT,1);
if (size(GT,2)~=size(HT,2))
    error('drsolve:t2cl:nargin','too few arguments')
end
a  = nroots1(m,xi);
b  = nroots1(n,eta);
GC = ftimes(GT,'A',xi);
HC = ftimes(HT,'A',eta);
