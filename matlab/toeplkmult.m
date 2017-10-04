function y = toeplkmult(G, B, x, ctrans, alg)
% y = toeplkmult(G, B, x)
% y = toeplkmult(G, B, x, ctrans)
% y = toeplkmult(G, B, x, ctrans, alg)
%
% Matrix-vector product Tx for a Toeplitz-like matrix T or T' with dense vectors x.
%
% Input:
%   G, B   -- Zp/Zm generators for T
%   x      -- dense vector(s)
%   ctrans -- (optional) compute T'*x instead of T*x, true / false (default)
%   alg    -- (optional) Algorithm choice
%              'full': Reconstruction of T + standard product (O(dn^2)).
%              'fft': Circulant based multiplication (O(d n log n), default).
%
% Output:
%   y      -- Dense vector(s) T*x (or T'*x)
%
% A concise description of the circulant based multiplication is
%
%   Bini CIMS paper TODO
%
% and it is also treated in Chapter 2 of
%
%   Pan book 2001 TODO.

if nargin < 5 || isempty(alg)
    alg = 'fft';
end

if nargin < 4 || isempty(ctrans)
    ctrans = false;
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
        if ~ctrans
            y = naive_alg(G,B,x);
        else
            y = naive_alg_ctrans(G, B, x);
        end
    case 'fft'
        % Basic blocking and avoidance of a few ifft calls
        if ~ctrans
            y = fftmult(G,B,x);
        else
            y = fftmult_ctrans(G,B,x);
        end
    case 'full'
        T = toeplkreconstruct(G,B);
        if ctrans
            y = T' * x;
        else
            y = T*x;
        end
    otherwise
        error('Invalid alg choice');
end

end

function y = naive_alg_ctrans(G, B, x)
% Naive implementation of the ctransposed multiplication.
% We use the fact that if
%
%   Zp * A  - A  * Zm = G * B',
%
% then A' satisfies the displacement equation
%
%   Zm * A' - A' * Zp = (Zm * B) * (Zp' * G)' =: GG * BB'.
%
% Hence we can use a similar circulant representation for A' as we did for
% A, viz.
%
%   -2 * T = sum_k C_m(gg_k) C_p(J * conj(bb_k))
%
% Where C_p(v) C_m(v) are f-circulant with first column v, and f=+1/-1.
%
% NOTE: This functions is a DEBUG function optimized for readability,
% and not for speed.
%
% See also: Pan (2001), Ch. 4.4

[n,r] = size(G);
nb = size(x,2);

GG = full(fcirculant(n,-1) * B); % the 'full' part affects only the 1x1 case.
BB = full(fcirculant(n,1)' * G); % idem

y = zeros(n,nb);

for k = 1:r
    tmp = circmult_p1(conj(BB(end:-1:1, k)), x);
    y = y + circmult_m1(GG(:,k), tmp);
end
y = -0.5 * y;
end

function y = naive_alg(G, B, x)
% Naive implementation, using the representation
%
%   2 * T = sum_k C_p(g_k) C_m(J*conj(b_k))
%
% Where C_p(v) C_m(v) are f-circulant with first column v, and f=+1/-1.
%
% NOTE: This functions is a DEBUG function optimized for readability,
% and not for speed.
%
% See also: Pan (2001), Ch. 4.4

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

% n-th unit root of f=-1
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
% See corresponding naive fft function for explanation.

realdata = isreal(G) && isreal(B) && isreal(x);

n = size(G, 1);
drank = size(G, 2);
nrhs = size(x,2);

% n-th unit root of f=-1
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

function y = fftmult_ctrans(G, B, x)
% See corresponding naive fft function for explanation.

realdata = isreal(G) && isreal(B) && isreal(x);

n = size(G, 1);
drank = size(G, 2);
nrhs = size(x,2);

% n-th unit root of f=-1
nrf = exp(1i * pi / n);

% Scaling matrix for the diagonalization of C_{-1}
d = transpose(nrf.^(0:n-1));
dinv = 1.0 ./ d;

% Transformation of generators from Zp/Zm -> Zm/Zp for A'
GG = full(fcirculant(n,-1) * B); % the 'full' part affects only the 1x1 case.
BB = full(fcirculant(n,1)' * G); % idem

% Transform the generators to FFT space, suiteable for circulant
% representation.
GG = fft(d .* GG);
BB = fft(conj(BB(end:-1:1, :)));

% Pull RHS vectors to fft space
x = fft(x);

% Accumulate product sum y in fftspace
y = zeros(n, nrhs);

for k=1:drank
    tmp = ifft(BB(:,k) .* x);
    if realdata
        % Gives a little bit of performance
        tmp = real(tmp);
    end
    y = y + GG(:,k) .* fft(d .* tmp);
end

y = -0.5 * dinv .* ifft(y);

if realdata
    y = real(y);
end


end

