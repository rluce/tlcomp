function [G, B] = toeplkgen(A, r)
% [G,B] = toeplkgen(A)
%
% Compute an approximate Z1/Zm1-generator for A.  If a rank r is supplied, the
% generator is truncated to that rank.
%
% CAUTION: expensive comutation, debug use only!
% TODO: We should do an ACA here instead of SVD, which would give us O(rn^2)
% run time.

n = size(A, 1);

if nargin >= 2 && ~isempty(r)
    if r > n
        error('tlcomp:InconsistentInput', ...
            'Requested rank exceeds matrix dimension');
    end
end

D = displace(A);
[U, S, V] = svd(D);

r_svd = rank(S);

if nargin < 2 || isempty(r)
    r = r_svd;
end
s = diag(S);

G = U(:,1:r) * diag(s(1:r));
B = V(:,1:r);


end
