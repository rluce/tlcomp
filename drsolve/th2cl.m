function [G,H,a,b] = th2cl(cT,rT,cH,rH)
%TH2CL	Converts a square Toeplitz+Hankel matrix to Cauchy-like.
%   [G,H,a,b] = T2CL(cT,rT,cH,rH) returns the quantities used in the
%   displacement equation
%
%       diag(a) * C - C * diag(b) = G * H'
%
%   of a Cauchy-like matrix C corresponding to the Toeplitz+Hankel
%   matrix M=T+H, where T=TOEPLITZ(cT,rT) and H=HANKEL(cH,rH).
%   The matrices C and M are such that C = U*M*V, where U=STIMES(n)
%   and V=CTIMES(n) are unitary matrices.
%
%   See also thl2cl, ctimes, stimes, toeplitz, hankel.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 19, 2010

if nargin<4, error('drsolve:th2cl:nargin','too few arguments'), end

cT = cT(:);
rT = flipud(rT(:)); % counterclockwise ordering in col vectors
cH = cH(:);
rH = rH(:);

n = size(cT,1);
if (n~=size(rT,1)) || (n~=size(cH,1)) || (n~=size(rH,1)) 
    error('drsolve:th2cl:nargin','wrong size')
end

if rT(n) ~= cT(1)
    warning('drsolve:th2cl:DiagonalConflict',['Toeplitz - First element ' ...
        'of input column does not match first element of input row. ' ...
        '\n         Column wins diagonal conflict.'])
end
rT(n) = cT(1);

if rH(1) ~= cH(n)
   warning('drsolve:th2cl:AntiDiagonalConflict',['Hankel - Last element ' ...
           'of input column does not match first element of input row. ' ...
           '\n         Column wins anti-diagonal conflict.'])
end
rH(1) = cH(n);

z = zeros(n-1,1);
G = [ ...
      cH-[0;cH(1:n-1)]+cT-[cT(2:n);0] ...
      [-1;z] ...
      [z;-1] ...
      rH-[rH(2:n);0]+rT-[0;rT(1:n-1)] ...
    ];

H = conj([ ...
      [-1;z] ...
      [0;cH(1:n-1)]+[flipud(rT(1:n-1));0] ...
      [rH(2:n);0]+[0;flipud(cT(2:n))] ...
      [z;-1]
    ]);

a = 2*cos((1:n).'*(pi/(n+1)));
b = 2*cos((0:n-1).'*(pi/n));

G = stimes(G);% ~dst
H = ctimes(H,'T');% dct
