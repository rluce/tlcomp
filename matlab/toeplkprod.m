function [Gp, Bp] = toeplkprod(G1, B1, G2, B2, alg)
% [Gs, Bs] = toeplkprod(G1, B1, G2, B2, alg)
%
% Compute generator for product Tp = T1 * T2 of two TL matrices.
%
% Input:
%   G1, B1 -- Generator for the TL matrix T1
%   G2, B2 -- Generator for the TL matrix T2
%   alg    -- Algorithmic choice:
%            'full': traditional matrix vector product with reconstructed
%                    matrix. Cost O((d1+d2)n^2).
%            'fft':  (default) FFT based matrix vector product.
%                    Cost O(d1 * d2 * n * log n).

% Generator formulas follow from the general Schur complement formula
% applied to the embedding
%
%   M = [
%         -I,  T2;
%         T1,   0;
%       ];
%
% This is spelled out directly in the 'full' code path.

if nargin < 5 || isempty(alg)
    alg = 'fft';
end

if isempty([G1, B1, G2, B2])
    Gp = [];
    Bp = [];
    return
end

switch alg
    case 'full'
        [Gp, Bp] = tlprod_full(G1,B1,G2,B2);
    case 'fft'
        [Gp, Bp] = tlprod_fft(G1,B1,G2,B2);
    otherwise
        error('Invalid choice for alg parameter');
end

end

function [Gp, Bp] = tlprod_fft(G1, B1, G2, B2)
n = size(G1, 1);

e1 = zeros(n,1);
e1(1) = 1;
en = zeros(n,1);
en(end) = 1;

tmp = toeplkmult(G1,B1,[e1, G2], false, 'fft');
Gp = [tmp, G1];
tmp = toeplkmult(G2, B2, [en, B1], true, 'fft');
Bp = [-2 * tmp(:,1), B2, tmp(:, 2:end)];
end

function [Gp, Bp] = tlprod_full(G1, B1, G2, B2)

n = size(G1, 1);

T1 = toeplkreconstruct(G1,B1);
T2 = toeplkreconstruct(G2,B2);

e1 = zeros(n,1);
e1(1) = 1;
en = zeros(n,1);
en(end) = 1;

Gp = [T1 * e1, T1 * G2, G1];
Bp = [-2 * T2' * en, B2, T2' * B1];

end
