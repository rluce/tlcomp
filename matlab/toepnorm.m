function val = toepnorm(c, r, p)
% val = toepnorm(c, r, p)  --  Compute matrix norm |T|_p of a Toeplitz
% matrix T.
%
% Input:
%   c, r -- first column and row of a Toeplitz matrix T
%   p    -- which norm to compute, possible choices:
%             1    : Matrix-1 norm (maximum column abs sum)
%             inf  : Matrix-inf norm (maximum row abs sum)
%             'fro': Frobenius norm
%
% Output:
%   val   -- |T|_p
%
% All of these norms are computed in O(n) ops and memory.
%
% The 2-norm cannot be computed as cheaply as the norms above, which is
% why it is not supported.  For this norm consider using |normest|.

% Homogenization
c = c(:);
r = r(:);

switch p
    case 1
        val = toep_norm_1(c,r);
    case inf
        val = toep_norm_1(conj(r), conj(c));
    case 'fro'
        val = toep_norm_fro(c,r);
    otherwise
        error('tlzstein:Unsupported', ...
            'Only matrix norms 1, inf and ''fro'' are supported');
end
end

function tnorm = toep_norm_fro(c,r)
n = length(c);
t = [r(end:-1:1); c(2:end)];
w = [1:n, (n-1):-1:1];
tnorm = w * abs(t).^2;
tnorm = sqrt(tnorm);
end

function tnorm = toep_norm_1(c,r)
% Compute 1-norm of a given Toeplitz matrix in O(n)

n = length(c);
colnorm = sum(abs(c));
tnorm = colnorm;

for j = 2:n
    colnorm = colnorm - abs(c(end-j+2)) + abs(r(j));
    tnorm = max(tnorm, colnorm);
end

end