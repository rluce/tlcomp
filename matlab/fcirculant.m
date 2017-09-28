function C = fcirculant(n, f)
% C = circulant(n, f)
%
% Returns the so-called f-circulant matrix of size n, i.e., for n=4
%   0 0 0 f
%   1 0 0 0
%   0 1 0 0
%   0 0 1 0
%
% By default f=0, so that C is the downshift matrix.

if nargin < 2 || isempty(f)
    f = 0.0;
end

if n <= 0
    % Convention: Empty f-circulant is empty+sparse
    C = sparse(0,0);
    return;
end

v = ones(n-1,1);
C = spdiags(v, -1, n, n);
C(1,n) = f;

end