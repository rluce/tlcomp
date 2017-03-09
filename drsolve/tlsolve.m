function [x rcond] = tlsolve(GT,HT,b,piv,xi,eta)
%TLSOLVE Solution of a Toeplitz-like linear system.
%   x = TLSOLVE(G,H,b) computes the solution to the linear system T*x=b,
%   where T is a square nonsingular Toeplitz-like matrix with
%   displacement equation
%
%      Z(1) * T - T * Z(-1) = G * H'.
%
%   b may contains multiple columns.
%   TLSOLVE converts the Toeplitz-like system to a Cauchy-like system,
%   and then solves it by CLSOLVE.
%
%   x = TLSOLVE(G,H,b,piv) calls CLSOLVE with the parameter piv, which
%   selects a pivoting technique. The default is piv=1 (partial
%   pivoting). See CLSOLVE for further details on pivoting.
%
%   [x,rcond] = TLSOLVE(...) also returns an estimate for the reciprocal
%   of the condition number given by TLSOLVE.
%
%   x = TLSOLVE(G,H,b,piv,xi,eta) represents the matrix T by the
%   displacement equation
%
%      Z(xi) * T - T * Z(eta) = GT * HT',           |xi|=|eta|=1,
%
%   where Z(xi) is the xi-cyclic forward shift matrix.
%
%   See also clsolve, tsolve.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

if nargin<3, error('drsolve:tlsolve:nargin','too few arguments'), end
if (nargin<4) || isempty(piv), piv = 1; end
if nargin<6, xi = 1; eta = -1; end

reale = isreal(GT) && isreal(HT) && isreal(xi) && isreal(eta) && isreal(b);

[GC,HC,tC,sC] = tl2cl(GT,HT,xi,eta);

bC        = ftimes(b,'A',xi);
[x rcond] = clsolve(GC,HC,tC,sC,bC,piv);
x         = ftimes(x,'N',eta);

if reale, x = real(x); end

%if rcond < eps
%	warning('drsolve:tlsolve:badConditioning',...
%    ['Matrix is close to singular or badly scaled.\n' ...
%    '         Results may be inaccurate. RCONDU = %e.'], rcond)
%end

