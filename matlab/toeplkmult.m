function y = toeplkmult(G, B, x, alg)
% y = toeplkmult(G, B, x)
%
% Multiplication of a Toeplitz-like matrix with a vector.

if nargin < 4 || isempty(alg)
    alg = 'full';
end

%y = tlmult_full(G, B, x);
%return;

switch alg
    case 'old'
        y = naive_alg(G,B,x);
    case 'new'
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

function y = naive_alg(G,B,x)
% Naive implementation, using the L_1 U_1 + ... + L_r U_r = T
% reconstruction formula.  Blocked only through the multiple RHS x.

[n,r] = size(G);
nb = size(x,2);

z = zeros(n,1);
y = zeros(n,nb);

for k = 1:r
    z(1) = conj(B(1,k));
    app_u = toepmult(z, conj(B(:,k)), x);
    z(1) = G(1,k);
    y = y + toepmult(G(:,k), z, app_u);
end

end
