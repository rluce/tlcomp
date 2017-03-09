function y = ctimes(x,s)
%CTIMES Cosine matrix product.
%   CTIMES multiplies the orthogonal cosine matrix C times a matrix X,
%   having at least two rows.
%
%   Y = CTIMES(X)         computes C*X;
%   Y = CTIMES(X,[])      as above;
%   Y = CTIMES(X,'N')     as above;
%   Y = CTIMES(X,'A')     computes C'*X (applies the adjoint of C);
%   Y = CTIMES(X,'T')     computes C.'*X (applies the transpose of C).
%
%   If the first parameter is a scalar, it must be a positive integer
%   and the output is a matrix, according to the second argument.
%  
%   Y = CTIMES(n)         returns C, which is n-by-n;
%   Y = CTIMES(n,[])      as above;
%   Y = CTIMES(n,'N')     as above;
%   Y = CTIMES(n,'A')     returns  C', i.e., the adjount of C;
%   Y = CTIMES(n,'T')     returns  C.', i.e., the transpose of C.
%
%   In Matlab notation C is shortly defined as C=idct(eye(n)).

%   See also ftimes, stimes, dct, idct.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

%==========================================================================
% parse input
%==========================================================================
% flag1: really assemble C [or C.']?
if size(x,1)>1
    n     = size(x,1);
    flag1 = 0;          % NO, I will apply [I]DCT of length n instead :)
elseif ~isscalar(x) || x~=ceil(x) || x<2
    error('ctimes:wrongX','x must be either n-by-m (and n>1), or an integer >1.');
else
    n     = x;
    flag1 = 1;          % YES, I will assemble a matrix :(
end
%--------------------------------------------------------------------------
% Traspose or Not?
if nargin<2 || isempty(s), s='N'; end
if s~='N' && s~='T' && s~='A'
    error('ctimes:wrongS','S must be ''T''ranspose [or ''A''] or ''N''ot');
end
%==========================================================================
% compute
%==========================================================================
if flag1
    %                                                       ASSEMBLE matrix
    if s=='N'
        y = idct(eye(n));%                                              C*I
    else
        % S='T' or 'A'
        y = dct(eye(n));%                                             C.'*I
    end
else
    %                                                          APPLY [i]dct
    %x=vector|matrix
    if s=='N'
        y = idct(full(x));%                                             C*x
    else
        % S='T' or 'A'
        y = dct(full(x));%                                            C.'*x
    end
end
