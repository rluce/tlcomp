function T = toeplkreconstruct(G,B)
% T = toeplktronstruct(G,B)
%  -- Reconstruct a Toeplitz-like matrix T from a generator (G,B).
%
% Input:
%   G, B  -- Generator of T w.r.t. Z1/Zm1-Stein displacement
%
% Output:
%   T     -- Reconstructed Toeplitz-like matrix
%
% Complexity O(d*n^2), where d is the displacement rank of T.


% Idea: Reconstruction by wrap-around cumulative sums and differences
n = size(G,1);
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