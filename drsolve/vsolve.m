function [x rcond] = vsolve(y,eta,b,piv)
%VSOLVE	Solution of a Vandermonde linear system.
%   x = VSOLVE(y,eta,b) computes the solution to the linear system
%   V*x=b where V = VANDER(y) is a square nonsingular matrix. The
%   parameter eta must be a complex number of modulus one, which is
%   different from y(i)^n, i=1,..,n. b may contain multiple columns.
%   VSOLVE converts the Vandermonde system to a Cauchy-like system,
%   and then solves it by CLSOLVE.
%
%   x = VSOLVE(y,eta,b,piv) calls CLSOLVE with the parameter piv,
%   which selects a pivoting technique. The default is piv=1 (partial
%   pivoting). See CLSOLVE for further details on pivoting.
%
%   [x,rcond] = VSOLVE(...) also returns an estimate for the
%   reciprocal of the condition number given by CLSOLVE.
%
%   See also clsolve, vlsolve.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

if nargin<3, error('drsolve:vsolve:nargin','too few arguments'), end
if nargin<4, piv = 1; end
y = y(:);
n = size(y,1);
if n~=size(b,1)
    error('drsolve:vsolve:size','wrong argument size.')
end

reale = isreal(y) && isreal(b);

[GC,HC,tC,sC] = v2cl(y,eta);

[x rcond] = clsolve(GC,HC,tC,sC,b,piv);
x         = ftimes(x,'N',eta);

if reale, x = real(x); end

%if rcond < eps
%	warning('drsolve:vsolve:badConditioning',...
%    ['Matrix is close to singular or badly scaled.\n' ...
%    '         Results may be inaccurate. RCONDU = %e.'], rcond)
%end

