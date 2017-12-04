function T = nice_spd_toeplitz(n)
% T = nice_spd_toeplitz(n)  --  SPD ToepMat with spectrum in [1e-1, 1e2].

T = ToepMat(randn(n,1));

opts.tol = 1e-5;
opts.isreal = true;
opts.issym = true;
opts.p = 40;
e = eigs(@(x) T*x, size(T,1), 2, 'BE', opts);
e_min = min(e);
e_max = max(e);

mu = (e_max - e_min) / (100 - 0.1);
sigma = 0.1 * mu - e_min;

T = T + sigma * toepeye(n);
T = T / mu;

end