function benchmark_tlmult(d, b)
% Compare performance of dense/fft multiplication

nmax = 2*4096;
dims = 128:128:nmax;

G = randn(nmax, d);
B = randn(nmax, d);
RHS = randn(nmax, b);

time_full = zeros(length(dims),1);
time_fft = zeros(length(dims),1);

for k = 1:length(dims)
    n = dims(k);
    tic;
    y_full = toeplkmult(G(1:n, :), B(1:n, :), RHS(1:n, :), 'full');
    time_full(k) = toc;
    
    tic;
    y_fft = toeplkmult(G(1:n, :), B(1:n, :), RHS(1:n, :), 'fft_naive');
    time_fft(k) = toc;
    
    fprintf('n=%4d  tfull=%2.2fs  tfft=%2.2fs  err=%.2e\n', ...
        n, time_full(k), time_fft(k), norm(y_full - y_fft)/norm(y_full));
end

plot(dims, [time_full, time_fft]);
legend({'full', 'fft'});
end