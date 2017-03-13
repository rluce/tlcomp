function [Gsylv, Bsylv] = stein2sylv(Gstein, Bstein)

n = size(Gstein, 1);

Zm1 = fcirculant(n, -1);

[c,r] = last_cross(Gstein, Bstein);
c = - circshift(c, 1, 1);
r = circshift(r, 1, 1);
r(1) = - r(1);

[Gadd, Badd] = stein_generator(c,r);
    
Gsylv = [-Gstein, Gadd];
Bsylv = [Zm1' * Bstein, Zm1' * Badd];

end

function [c,r] = last_cross(G,B)
% Compute last row and column of a Toeplitz-like matrix

n = size(G,1);
GB = G*B';
c = zeros(n,1);
r = zeros(n,1);

% Summing subdiagonals
for k = -(n-1):0
    r(k + n) = sum(diag(GB,k));
end

% Summing superdiagonals
for k = 0:(n-1)
    c(n-k) = sum(diag(GB,k));
end

end