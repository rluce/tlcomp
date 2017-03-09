function [c, r, T] = random_toeplitz(m, n)
% [c, r, T] = random_toeplitz(m, n)
%
% Generate a random, complex Toeplitz matrix.

z = randn(m+n-1,1) + 1i * randn(m+n-1,1);
c = z(m:-1:1);
r = z(m:end).';
T = toeplitz(c,r);

end