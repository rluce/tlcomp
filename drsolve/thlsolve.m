function [x rcond] = thlsolve(G,H,b,piv)
%THLSOLVE Solution of a Toeplitz+Hankel-like linear system.
%   x = THLSOLVE(G,H,b) computes the solution to the linear system
%   Mx=b, where M is a square nonsingular Toeplitz+Hankel-like matrix
%   with displacement equation
%
%      Y_0 * M - M * Y_1 = G * H',
%
%   where Y_0 = tridiag(1,0,1) and Y_1 equals Y_0 but has 1 in the
%   'top left' and 'bottom right' corners. b may contain multiple
%   columns.
%   THLSOLVE converts the Toeplitz+Hankel-like system to a Cauchy-like
%   system, and then solves it by CLSOLVE.
%
%   x = THLSOLVE(G,H,b,piv) calls CLSOLVE with the parameter piv,
%   which selects a pivoting technique. The default is piv=1 (partial
%   pivoting). See CLSOLVE for further details on pivoting.
%
%   [x,rcond] = THLSOLVE(...) also returns an estimate for the reciprocal 
%   of the condition number given by CLSOLVE.
%
%   See also clsolve, thsolve, tsolve, tlsolve.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

if nargin<3, error('drsolve:thlsolve:nargin','too few arguments'), end
if nargin<4, piv = 1; end

n = size(G,1);
if (n~=size(H,1)) || (size(G,2)~=size(H,2))|| (n~=size(b,1)) 
    error('drsolve:thlsolve:size','check input arguments.')
end

reale = isreal(G) && isreal(H) && isreal(b);

[GC,HC,tC,sC] = thl2cl(G,H);
bC            = stimes(b);% ~dst
[x rcond]     = clsolve(GC,HC,tC,sC,bC,piv);
x             = ctimes(x,'N');% idct

if reale, x = real(x); end

%if rcond < eps
%	warning('drsolve:thlsolve:badConditioning',...
%    ['Matrix is close to singular or badly scaled.\n' ...
%    '         Results may be inaccurate. RCONDU = %e.'], rcond)
%end

