function [x rcond] = thsolve(cT,rT,cH,rH,b,piv)
%THSOLVE Solution of a Toeplitz+Hankel linear system.
%   x = THSOLVE(cT,rT,cH,rH,b) computes the solution to the linear
%   system (T+H)*x=b, where T = TOEPLITZ(cT,rT) and H = HANKEL(cH,rH)
%   are square matrices; b may contain multiple columns.
%   THSOLVE converts the Toeplitz-plus-Hankel system to a Cauchy-like
%   system, and then solves it by CLSOLVE.
%
%   x = THSOLVE(cT,rT,cH,rH,b,piv) calls CLSOLVE with the parameter
%   piv, which selects a pivoting technique. The default is piv=1
%   (partial pivoting). See CLSOLVE for further details on pivoting.
%
%   [x,rcond] = THSOLVE(...) also returns an estimate for the reciprocal 
%   of the condition number given by CLSOLVE.
%
%   See also clsolve, thlsolve, tsolve, tlsolve.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

if nargin<5, error('drsolve:thsolve:nargin','too few arguments'), end
if nargin<6, piv = 1; end
cT = cT(:);
rT = rT(:);
cH = cH(:);
rH = rH(:);
n = size(cT,1);
if (n~=size(rT,1)) || (n~=size(rH,1))|| (n~=size(cH,1)) || (n~=size(b,1)) 
    error('drsolve:thsolve:size','check input arguments.')
end

reale = isreal(cT) && isreal(rT) && isreal(cH) && isreal(rH) && isreal(b);

[GC,HC,tC,sC] = th2cl(cT,rT,cH,rH);
bC            = stimes(b);% ~dst
[x rcond]     = clsolve(GC,HC,tC,sC,bC,piv);
x             = ctimes(x,'N');% idct

if reale, x = real(x); end

%if rcond < eps
%	warning('drsolve:thsolve:badConditioning',...
%    ['Matrix is close to singular or badly scaled.\n' ...
%    '         Results may be inaccurate. RCONDU = %e.'], rcond)
%end

