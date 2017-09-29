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
    alg = 'full';
end

%y = tlmult_full(G, B, x);
%return;

switch alg
    case 'fft_naive'
        y = naive_alg(G,B,x);
    case 'fft'
        y = lu_circulant_expand(G,B,x);
    case 'full'
        y = tlmult_full(G, B, x);
    otherwise
        error('Invalid alg choice');
end

end

function y = tlmult_full(G, B, x)
T = toeplkreconstruct(G,B);
y = T*x;
end

function y = lu_circulant_expand(G,B,x)
[n,r] = size(G);
nb = size(x,2);

N = 2*n - 1;

% FFT of X
Fx = fft(x, N, 1);

% FFT of U
FU = fft([conj(B(1,:)); zeros(n-1, r); conj(B(end:-1:2,:))], N);

% FFT of L
FL = fft(G, N, 1);

% Accumulation of the (transformed) terms L_k * U_k * X
acc = zeros(N, nb);

for k=1:r
    % TODO Is there a way to avoid the fft/ifft transforms inside this
    % loop?  Maybe there exists a clever circulant embdedding for L*U that
    % does the trick?
    
    % In later Matlab versions, something goes wrong for singletons here,
    % i.e. if size(G,1) == 1 so that T is a scalar.  Unclear what the best
    % fix is. TODO.
    tmp = repmat(FU(:,k), 1, nb) .* Fx;
    tmp = ifft(tmp);
    tmp = tmp(1:n,:);
    tmp = fft(tmp,N,1);
    tmp = repmat(FL(:,k), 1, nb) .* tmp;
    acc = acc + tmp;
end

% IFFT tranform of result
Y = ifft(acc);
y = Y(1:n,:);

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

