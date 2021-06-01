function d = toeplkdet(G, B)
% d = toeplkdet(G, B) -- Determinant of a TL matrix.
%
% Input:
%   G, B -- Generator of the TL matrix A
%
% Output:
%   d    -- det(A)
%
% NOTE: The computation is based on the LU factorization of A after FFT to
% a Cauchy matrix.  The limitations for the usual det function apply here as
% well.

% The factorization we compute has the form
%
%   A = F' * P * L * U * F * D,
%
% where F is the normalized FFT, (P,L,U) the LU factorization of the TL
% matrix in FFT space, and D the diagnal matrix with the n unit roots of
% 1.0.

n = size(G,1);

if n == 1
    % Corner cases would mees up syntax for FFT, do it directly.
    d = toeplkreconstruct(G,B);
    return;
end

[GC,HC,a,b] = tl2cl(G, B, 1, -1);
d = gkolu(a, b, GC, HC');

% Account for the first n unit roots
root_diag = exp( ((0:n-1) * 1i * pi) / n );
d = prod(root_diag) * d;
if isreal(G) && isreal(B)
    d = real(d);
end

end

function d = gkolu(s, t, G, B)
% Determinant of the Cauchy matrix w.r.t. even/odd 2nth unit roots.

n = size(G,1);
p = 1:n;
diag_U = zeros(n,1);

% Track the signum of the permutation matrix
permsign = 1.0;

for k=1:n
    % Get first column of Schur complement, chose absmax row
    l(k:n) = (G(k:n,:)*B(:,k))./(s(k:n)-t(k));
    [abspiv, pivotpos] = max(abs(l(k:n)));
    
    if abspiv == 0.0
        diag_U(k) = 0.0;
        break;
    end

    if pivotpos ~= 1
        permsign = -1 * permsign;
    end

    % Permute G generator and s data, track permutation
    pivotpos = pivotpos + k - 1;
    p([k, pivotpos]) = p([pivotpos, k]);
    G([k, pivotpos],:) = G([pivotpos, k],:);
    s([k, pivotpos]) = s([pivotpos, k]);

    % Compute next col/row of L/U
    u = (G(k,:)*B(:,k:n))./ transpose(s(k)-t(k:n));
    l = (G(k+1:n,:)*(B(:,k)/u(1))) ./ (s(k+1:n)-t(k));
    diag_U(k) = u(1);
    
    % Compute generator for Schur complement
    G(k+1:n,:) = G(k+1:n,:) - l * G(k,:);
    B(:,k+1:n) = B(:,k+1:n) - (B(:,k)/u(1,1)) * u(2:end);
end

d = permsign * prod(diag_U);
end
