function benchmark_tlmult(d, b)
% Compare performance of dense/fft multiplication

nmax = 2*4096;
dims = 128:128:nmax;

G = randn(nmax, d) + 1i * randn(nmax, d);
B = randn(nmax, d) + 1i * randn(nmax, d);
RHS = randn(nmax, b) +1i * randn(nmax, b);

time_full = zeros(length(dims),1);
time_fft = zeros(length(dims),1);
time_fft2 = zeros(length(dims),1);

for k = 1:length(dims)
    n = dims(k);
    tic;
    y_full = toeplkmult(G(1:n, :), B(1:n, :), RHS(1:n, :), 'full');
    time_full(k) = toc;
    
    tic;
    y_fft = toeplkmult(G(1:n, :), B(1:n, :), RHS(1:n, :), 'fft_naive');
    time_fft(k) = toc;

    tic;
    y_fft2 = toeplkmult(G(1:n, :), B(1:n, :), RHS(1:n, :), 'fft');
    time_fft2(k) = toc;

    fprintf('n=%4d  tfull=%2.2fs  tfft=%2.2fs  tfft2=%2.2fs  err1=%.2e  err2=%.2e\n', ...
        n, time_full(k), time_fft(k), time_fft2(k), ...
        norm(y_full - y_fft)/norm(y_full), norm(y_full - y_fft2)/norm(y_full));
end

plot(dims, [time_full, time_fft, time_fft2]);
legend({'full', 'fft', 'fft2'});
end