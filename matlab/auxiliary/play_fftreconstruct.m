function pgnwiybeg

n = 8;
f = -1;

% First n-th root of  f (f is on the unit circle
nrf = exp(1i * angle(f) / n);
v = (1:n)';
v = randn(n,1);
d = transpose(nrf.^(0:n-1));

Z = fcirc(f, v)
Omega = fft(eye(n));

ZZ = 1/n * diag(1.0./d) * Omega' * diag(Omega * diag(d) * v) * Omega * diag(d)

norm(Z - ZZ)

x = randn(n,1);
y = Z * x;
yy = real((1.0./d) .* ifft( fft(d.*v) .* fft(d.*x) ));
norm(y - yy)

% Basic reconstruction formula test
n = 9;
G = orth(randn(n,1));
B = orth(randn(n,1));

T_ref = toeplkreconstruct(G,B);

% Sum up products of circulants
T = zeros(n);
for k = 1:size(B,2)
    T = T + fcirc(1, G(:,k)) * fcirc(-1, conj(B(end:-1:1,k)));
end
T = 0.5 * T;

norm(T_ref - T)
end

