function [Gs, Bs] = toeplksquare(G, B, alg)
% [Gs, Bs] = toeplksquare(G, B, alg)
%
% Compute generator for squared TL matrix.
%
% Possible choices for paramter 'alg':
%   'full' -- traditional matrix vector product with reconstructed matrix
%   'fft'  -- FFT based matrix vector product

if nargin < 3 || isempty(alg)
    alg = 'fft';
end

switch alg
    case 'full'
        [Gs, Bs] = tlsquare_full(G,B);
    case 'fft'
        [Gs, Bs] = tlsquare_fft(G,B);
    otherwise
        error('Invalid choice for alg parameter');
end

end

function [Gs, Bs] = tlsquare_fft(G,B)

n = size(G,1);
e = zeros(n,1);
e(1) = 1;
Gs = [ toeplkmult(G, B, [e, G], 'fft'), G];
e(1) = 0;
e(n) = 1;
tmp = toeplkmult(conj(B), conj(G), [e, B], 'fft');
Bs = [-2 * tmp(:,1), B, tmp(:, 2:end)];
end

function [Gs, Bs] = tlsquare_full(G, B)
T = toeplkreconstruct(G,B);
Gs = [ T(:,1), T*G, G ];
Bs = [ -2*T(end,:)', B, T' * B ];
end
