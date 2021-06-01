function quicksigntest


n = 1000;

t = randn(n,1);

T = ToepMat(t,t);
tic;
S = signm(full(T));
fprintf('Full signm(T): %.1fs\n', toc);
maxiter = 10;
tim = 0.0;
for k = 1:maxiter
    tic;
    T = 0.5 * (T + inv(T));
    tim = tim + toc;
    fprintf('Iteration %2d, drank=%2d (bound: %4d) error=%.2e  tim=%5.1fs\n', ...
        k, drank(T), 2^(k+1), norm(full(T) - S, 'fro') / norm(S, 'fro'), tim);

end

end

