function [G,H,xi,eta] = t2tl(c,r,phi,psi)
%T2TL	Converts a Toeplitz matrix to Toeplitz-like.
%   [G,H,xi,eta] = T2TL(c,r) returns the generators of the matrix
%   T=TOEPLITZ(c,r) with respect to the displacement equation
%
%      Z(xi) * T - T * Z(eta)  =  G * H'.
%
%   The sizes are:
%      matrix  |  T      Z(xi)  Z(eta)  G      H      xi     eta
%      size    |  [m,n]  [m,m]  [n,n]   [m,2]  [n,2]  [1x1]  [1x1]
%
%   Here Z(xi) is the square xi-cyclic forward shift matrix defined as
%   Z(xi) = [O,xi;eye(m-1),0], where 0 are null blocks.
%
%   [...] = T2CL(c,r,xi,eta) forces a choice for xi, eta. The default
%   is xi=1 and eta=exp(1i*pi*gcd(m,n)/m). The two scalars xi and eta
%   must be unitary, i.e., abs(xi)=abs(eta)=1.
%
%   See also t2cl, tl2cl, toeplitz.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

% T is m-x-n
% build t = [ T_{1,n} ... T_{1,1} ... T_{m,1} ].' ...
if nargin<2, error('drsolve:t2tl:nargin','too few arguments'), end
c = c(:);
r = r(:);
m = size(c,1);
n = size(r,1);
if r(1) ~= c(1)
    warning('drsolve:t2tl:DiagonalConflict',['First element of ' ...
        'input column does not match first element of input row. ' ...
        '\n         Column wins diagonal conflict.'])
end
t = [r(n:-1:2) ; c];                 % build vector "t" of user data

% define xi eta
if nargin==2
    xi = 1;
    if mod(n,m), eta = exp(complex(0,pi*(gcd(m,n)/m))); else eta = -1; end
    %
elseif nargin~=4 || abs(abs(phi)-1)>eps || abs(abs(psi)-1)>eps
    error('drsolve:t2tl:nargin','check input arguments')
else
    xi=phi; eta=psi; % set up for LHS
end
%if nargout<4 && nargin<4, warning('drsolve:t2tl:output','not injective'), end

% define G H
% use +n wrt formula :)
G = [ ... 
               -eta*t(n)           1;
       t(1:m-1)-eta*t(n+1:n+m-1)   zeros(m-1,1)
    ];
H = [ ...
       zeros(n-1,1)  conj(xi*t(n+m-1:-1:m+1)-t(n-1:-1:1));
       1             conj(xi*t(m))
    ];
