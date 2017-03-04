function [G,B] = compute_generator(A, r)
% [G,B] = compute_generator(A)
%
% Compute some Z-generator for A.  If r is supplied, the generator is
% truncated to that rank.
%
% CAUTION: expensive comutation, debug use only!

n = size(A, 1);

if nargin >= 2 && ~isempty(r)
    if r > n
        error('tlzstein:InconsistentInput', ...
            'Requested rank exceeds matrix dimension');
    end
end

Z = downshift(n);
D = A - Z * A * Z';
[U, S, V] = svd(D);

r_svd = rank(S);

if nargin < 2 || isempty(r)
    r = r_svd;
end
s = diag(S);

G = U(:,1:r) * diag(s(1:r));
B = V(:,1:r);


end
