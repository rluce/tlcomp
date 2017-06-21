function nveoiuybaawe

maxiter = 10;

r = randn(2500,1);
T = toeplitz(r);
e_min = eigs(T, 1, 'lm');
r(1) = r(1) + 0.01 + abs(e_min);


T = toeplitz(r);
e_min = eigs(T,1,'sm');
e_max = eigs(T,1,'lm');

fprintf('min/max eigenvalue: %.2e/%.2e\n', e_min, e_max);
tic;
S = sqrtm(T);
tim = toc;
normS = norm(S, 'fro');
fprintf('sqrtm took %.2fs, norm is %.2e\n', tim, normS);

TL = TLMat(r,r);
X = TL;
fprintf('%4s  %4s  %11s\n', 'Iter', 'dprk', 'relerr')
fprintf('%4d  %4d  %11.2e\n', 0, drank(X), norm(X - S, 'fro')/normS);

Xmu = X;

mu = zeros(1,maxiter);
kappa = sqrt(e_min/e_max);
mu(1) = (e_min * e_max)^(-1/4);
mu(2) = sqrt( 2 * sqrt(kappa) / (1 + kappa) );
for k = 3:maxiter
    mu(k) = sqrt( 2 * mu(k-1) / (1 + mu(k-1)^2) );
end

for k = 1:maxiter
    X = .5 * (X + TL * inv(X)); %#ok<*MINV>
    Xmu = .5* ( mu(k) * Xmu + (1.0/mu(k)) * TL * inv(Xmu));
    fprintf('%4d  %4d  %11.2e  %.2e  %4d  %11.2e\n', k, ...
        drank(X), norm(X - S, 'fro')/normS, mu(k), drank(Xmu), norm(Xmu - S, 'fro')/normS);
end

end