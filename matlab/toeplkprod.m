function [Gp, Bp] = toeplkprod(G1, B1, G2, B2, alg)
% [Gs, Bs] = toeplkprod(G1, B1, G2, B2, alg)
%
% Compute generator for product of TL matrices.
%
% Possible choices for paramter 'alg':
%   'full' -- traditional matrix vector product with reconstructed matrix
%   'fft'  -- FFT based matrix vector product

if nargin < 5 || isempty(alg)
    alg = 'full';
end

switch alg
    case 'full'
        [Gp, Bp] = tlsquare_full(G1,B1,G2,B2);
    case 'fft'
        [Gp, Bp] = tlsquare_fft(G1,B1,G2,B2);
    otherwise
        error('Invalid choice for alg parameter');
end

end

function [Gp, Bp] = tlsquare_fft(G1, B1, G2, B2)



n = size(G1,1);

onesvec = ones(n,1);

Gp_part1 = vapply(G2, 'inv');
Gp_part1 = toeplkmult(G1, B1, Gp_part1);
Gp_part1 = vapply(Gp_part1);

Gp_part2 = toeplkmult(G1, B1, -onesvec);
Gp_part2 = vapply(Gp_part2);

Gp = [Gp_part1, G1, -Gp_part2];

Bp_part1 = vapply(B1, 'inv');
Bp_part1 = toeplkmult(B2, G2, Bp_part1);
Bp_part1 = vapply(Bp_part1);

Bp_part2 = toeplkmult(B2, G2, -onesvec);
Bp_part2 = vapply(Bp_part2);

Bp = [B2, Bp_part1, Bp_part2];
end

function [Gp, Bp] = tlsquare_full(G1, B1, G2, B2)
T1 = stein_reconstruction(G1,B1);
T2 = stein_reconstruction(G2,B2);

Gp_part1 = vapply(G2, 'inv');
Gp_part1 = T1*Gp_part1;
Gp_part1 = vapply(Gp_part1);

Gp_part2 = - sum(T1,2);
Gp_part2 = vapply(Gp_part2);

Gp = [Gp_part1, G1, -Gp_part2];

Bp_part1 = vapply(B1, 'inv');
Bp_part1 = T2' * Bp_part1;
Bp_part1 = vapply(Bp_part1);

Bp_part2 = - sum(T2,1)';
Bp_part2 = vapply(Bp_part2);

Bp = [B2, Bp_part1, Bp_part2];
end
