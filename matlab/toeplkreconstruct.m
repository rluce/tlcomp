function T = toeplkreconstruct(G, B)
% T = toeplktronstruct(G, B)
%  -- Reconstruct a Toeplitz-like matrix T from a generator (G,B).
%
% Input:
%   G, B  -- Generator of T w.r.t. Z1/Zm1 Sylvester displacement
%
% Output:
%   T     -- Reconstructed, full Toeplitz-like matrix
%
% Complexity O(d*n^2), where d is the displacement rank of T.

% TODO: This function is used for toeplkmult, where not the matrix T is
% needed as a whole, but only for the purpose of computing T*x. In this
% setting the reconstruction technique used here
% could be optimized for memory consumption.  As of now the full
% outer product T=G*B' is computed first, and then processed in a second
% step.  But there T is only accessed only column-by-column, so we could
% generate one (or a small batch of) column, process it, use it and forget
% it.  Consequently we would use only O(n) memory instead of O(n^2).

% TODO: It's worth some time to figure out how to reconstruct sections of
% the matrix, in particular bands and submatrices.  Some of these cases can
% certainly be better implemented than constructing the whole matrix first.

n = size(G, 1);

if size(B, 2) ~= size(G, 2)
    % Generator matrices are not compatible
    error('tlcomp:InconsistentInput', ...
            'G and B must have the same number of columns');
end

% Only implemented for square, could possibly be lifted.
assert(size(B, 1) == n);

if n == 0
    T = [];
    return;
end

% General idea: Reconstruction by wrap-around cumulative sums and differences

T = G*B';

% Compute diagonal-wrap-around sums (Use sym D 4x4 to understand the result)
s = T(:,1);
for k = 2:n
    tmpval = s(end) + T(1,k);
    s(2:end) = s(1:end-1) + T(2:end,k);
    s(1) = tmpval;
end

% This is the first column of T
newvals = .5 * s;

% Remaining columns arise from first by
% cumulative-wrap-around-diagonal-sums
% Use sym D 4x4 to understand the result. HINT: you need to comment out the next
% line.
%tmpvec = zeros(n,1);
for k = 1:n
    % tmpvec will hold the updated values for column k+1
    tmpvec(1,1) = -T(1,k) + newvals(end); % Wrap around value
    tmpvec(2:n,1) = -T(2:end, k) + newvals(1:end-1); % No wrap in rest of col
    
    % Column k no longer needed, now it gets its final value
    T(:, k) = newvals;
    
    % Iteration upkeep: Carry over values for column k+1
    newvals = tmpvec;
end

end
