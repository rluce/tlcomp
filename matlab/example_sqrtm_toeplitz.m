function example_sqrtm_toeplitz
r = zeros(1000,1);
r(1) = 3;
r(2) = -1;
t_tot = tic;
S = sqrtm_toeplitz(r);
fprintf('Total time: %.2fs\n', toc(t_tot));

if length(r) <= 1000
    % Don't do it for large matrices..
    S_true = sqrtm(toeplitz(r));
    fprintf('Relative dist to Matlab sqrtm: %8.2e\n', ...
        norm(S - S_true, 'fro') / norm(S_true, 'fro'));
end

end