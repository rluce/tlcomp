function [x rcond] = vlsolve(y,eta,G,H,b,piv)
%VLSOLVE Solution of a Vandermonde-like linear system.
%   x = VLSOLVE(y,eta,G,H,b) computes the solution to the linear
%   system V*x=b where V is a square nonsingular Vandermonde-like
%   matrix with displacement equation
%
%      diag(y) * V - V * Z(eta) = G * H'.
%
%   The parameter eta must be a complex number of modulus one, which
%   is different from y(i)^n, i=1,..,n. b may contain multiple
%   columns.
%   VLSOLVE converts the Vandermonde system to a Cauchy-like system,
%   and then solves it with CLSOLVE.
%
%   x = VLSOLVE(y,eta,G,H,b,piv) calls CLSOLVE with the parameter piv,
%   which selects a pivoting technique. The default is piv=1 (partial
%   pivoting). See CLSOLVE for further details on pivoting.
%
%   [x,rcond] = VLSOLVE(...) also returns an estimate for the
%   reciprocal of the condition number given by CLSOLVE.
%
%   See also clsolve, vsolve.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

if nargin<5, error('drsolve:vlsolve:nargin','too few arguments'), end
if nargin<6, piv = 1; end

y = y(:);
n = size(y,1);
if n~=size(G,1) || n~=size(G,1) || n~=size(b,1) || size(G,2)~=size(H,2)
    error('drsolve:vlsolve:size','wrong argument size.')
end

reale = isreal(y) && isreal(b) && isreal(G) && isreal(H) && isreal(eta);

[GC,HC,tC,sC] = vl2cl(y,eta,G,H);

[x rcond] = clsolve(GC,HC,tC,sC,b,piv);
x         = ftimes(x,'N',eta);

if reale, x = real(x); end

%if rcond < eps
%	warning('drsolve:vlsolve:badConditioning',...
%    ['Matrix is close to singular or badly scaled.\n' ...
%    '         Results may be inaccurate. RCONDU = %e.'], rcond)
%end

