function y = stimes(x)
%STIMES Sine matrix product.
%   STIMES multiplies the symmetric orthogonal sine matrix S times a
%   matrix X, having at least two rows.
%
%   Y = STIMES(X)       computes S*X.
%
%   If the parameter is a scalar, it must be a positive integer and
%   the output is the matrix S;
%  
%   S = STIMES(n)       returns S, which is n-by-n.
%
%   In Matlab notation S is shortly defined as
%   S=sqrt(2/(n+1))*dst(eye(n)). 
%
%   See also ctimes, ftimes, dst.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Dec 11, 2009

%==========================================================================
% parse input
%==========================================================================
% flag1: build a matrix?
if size(x,1)>1
    n     = size(x,1);
    flag1 = 0;          % NO, I will apply DST of length n instead
elseif ~isscalar(x) || x~=ceil(x) || x<2
    error('stimed:wrongX','x must be either n-by-m (and n>1), or an integer >1.');
else
    n     = x;
    flag1 = 1;          % YES, I will assemble am matrix
end
%==========================================================================
% compute
%==========================================================================
if flag1
        y = sqrt(2/(n+1))*dst(eye(n));%                                 S*I
else
        y = sqrt(2/(n+1))*dst(full(x));%                                S*x
end
