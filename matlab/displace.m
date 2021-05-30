function D = displace(A)
% D = displace(A) -- compute displacement of a given matrix A.
%
% Input:
%   A -- input matrix
%
% Output:
%   D -- displacement of A, i.e., D = Z1 * A - A * Zm1.  D is full.
%
% Complexity: O(n^2)

[m, n] = size(A);
Z1 = fcirculant(m, 1);
Zm1 = fcirculant(n, -1);

% Sparsity important for performance.
assert(issparse(Z1) && issparse(Zm1));
D = Z1 * A - A * Zm1;

% Affects only m=1 or n=1, as in this case D is sparse, so we convert.
D = full(D);
end