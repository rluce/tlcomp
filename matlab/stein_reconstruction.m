function T = stein_reconstruction(G,B)
% T = stein_reconstruction(G,B)
%  -- Reconstruct a Toeplitz-like matrix T from a generator (G,B).
%
% Input:
%   G, B  -- Generator of T w.r.t. Z-Stein displacement
%
% Output:
%   T     -- Reconstructed Toeplitz-like matrix
%
% Complexity O(d*n^2), where d is the displacement rank of T.

n = size(G,1);

% Idea: Recover T by forming cumulative sums along all diagonals of G*B'.
T = G*B';

% First half of diagonals
for j = 1:(n-1)
    T((j+1):end, j+1) = T((j+1):end, j+1) + T(j:(end-1), j);
end

% Second half of diagonals.  Work on transpose to avoid strided memory
% access.
T = T';
for j = 1:(n-1)
    T((j+2):end, j+1) = T((j+2):end, j+1) + T((j+1):(end-1), j);
end
T = T';

end