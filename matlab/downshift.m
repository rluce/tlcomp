function C = downshift(n)
% C = downshift(n)
%
% Returns the downshift matrix of size n, i.e., for n=4
%   0 0 0 0
%   1 0 0 0
%   0 1 0 0
%   0 0 1 0


v = ones(n-1,1);
C = spdiags(v, -1, n, n);

end