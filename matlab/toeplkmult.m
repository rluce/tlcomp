function y = toeplkmult(G, B, x, alg)
% y = toeplkmult(G, B, x)
% y = toeplkmult(G, B, x, alg)
%
% Matrix-vector product Tx for a Toeplitz-like matrix T with dense vectors x.
%
% Input:
%   G, B -- Z1/Zm1 generators for T
%   x    -- dense vector(s) 
%   alg  -- (optional) Algorithm choice
%             'full': Reconstruction of T + standard product (O(dn^2)).
%             'fft': Circulant based multiplication (O(d n log n))
%
% Output:
%   y    -- Dense vector(s) T*x
%
% An concise description of the circulant based multiplication is
%
%   Bini CIMS paper TODO
%
% and it is also treated in Chapter 2 of
%
%   Pan book 2001 TODO.

if nargin < 4 || isempty(alg)
    alg = 'fft';
end

% Corner case: T is 1x1 and the generator is non-minimal.  This confuses
% the way the dimension in the FFT is inferred.  Instead of obfuscating the
% syntax in the blocked FFT algorithm, we redirect this corner case.
if size(B, 1) == 1 && size(B, 2) > 1 && strcmp(alg, 'fft')
    alg = 'fft_naive';
end

switch alg
    case 'fft_naive'
        % Readable implementation of the circulant expansion
        y = naive_alg(G,B,x);
    case 'fft'
        % Basic blocking and avoidance of a few ifft calls
        y = fftmult(G,B,x);
    case 'full'
        T = toeplkreconstruct(G,B);
        y = T*x;
    otherwise
        error('Invalid alg choice');
end

end

function y = naive_alg(G, B, x)
% Naive implementation, using the representation
%
%   2 * T = sum_k C_1(g_k) C_m1(b_k)
%
% Where C_1(v) C_m1(v) are f-circulant with first column v, and f=1/-1.
%
% NOTE: This functions is a DEBUG function optimized for readability,
% and not for speed.

[n,r] = size(G);
nb = size(x,2);

y = zeros(n,nb);

for k = 1:r
    tmp = circmult_m1(conj(B(end:-1:1, k)), x);
    y = y + circmult_p1(G(:,k), tmp);

end
y = .5 * y;
end

function y = circmult_p1(v,x)
% Matrix vector multiplication with C_1(v)
realdata = isreal(v) && isreal(x);

% Pull data to fft space
v = fft(v);
x = fft(x);

% Multiply in fft space
y = v .* x;

% Pull back to original coordinates
y = ifft(y);

if realdata
    y = real(y);
end

end

function y = circmult_m1(v,x)
% Matrix vector multiplication with C_{-1}(v)

realdata = isreal(v) && isreal(x);

n = length(v);

% n-th unit root of -1
nrf = exp(1i * pi / n);

% Scaling matrix for the diagonalization of C_{-1}
d = transpose(nrf.^(0:n-1));
dinv = 1.0 ./ d;

% Pull data to fft space
v = fft(d.*v);
x = fft(d.*x);

% Multiply in fft space
y = v .* x;

% Pull back to original coordinates
y = dinv .* ifft(y);

if realdata
    y = real(y);
end

end

function y = fftmult(G, B, x)

realdata = isreal(G) && isreal(B) && isreal(x);

n = size(G, 1);
drank = size(G, 2);
nrhs = size(x,2);

% n-th unit root of -1
nrf = exp(1i * pi / n);

% Scaling matrix for the diagonalization of C_{-1}
d = transpose(nrf.^(0:n-1));
dinv = 1.0 ./ d;

% Transform the generators to FFT space, suiteable for circulant
% representation.
G = fft(G);
B = fft(d .* conj(B(end:-1:1, :)));

% Pull scaled RHS vectors to fft space
x = fft(d .* x);

y = zeros(n, nrhs);

for k=1:drank
    tmp = dinv .* ifft(B(:,k) .* x);
    if realdata
        % Gives a little bit of performance
        tmp = real(tmp);
    end
    y = y + G(:,k) .* fft(tmp);
end

y = 0.5 * ifft(y);

if realdata
    y = real(y);
end

end

