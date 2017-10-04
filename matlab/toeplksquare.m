function [Gs, Bs] = toeplksquare(G, B, alg)
% [Gs, Bs] = toeplksquare(G, B, alg)
%
% Compute generator for squared TL matrix.
%
% Possible choices for paramter 'alg':
%   'full' -- traditional matrix vector product with reconstructed matrix
%   'fft'  -- FFT based matrix vector product

% Idea: Use the augmentation
%
%   M = [ -I, A;
%          A, 0 ]
%
% so that the Schur complement of (1,1) in M is A^2.
%
% From the generic generator formula on obtains that the generator is
% given through
%
%   Gs = [ T*e1,         T*G, G ]
%   Bs = [ -2 * T' * en, B,   T' * B]
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
Gs = [ toeplkmult(G, B, [e, G], false, 'fft'), G];
e(1) = 0;
e(n) = 1;
tmp = toeplkmult(G, B, [e, B], true, 'fft');
Bs = [-2 * tmp(:,1), B, tmp(:, 2:end)];
end

function [Gs, Bs] = tlsquare_full(G, B)
T = toeplkreconstruct(G,B);
Gs = [ T(:,1), T*G, G ];
Bs = [ -2*T(end,:)', B, T' * B ];
end
