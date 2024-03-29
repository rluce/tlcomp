function x = toepsolve(c, r, b)
% x = toepsolve(c, r, b)
%
% Solve a Toeplitz system of equations in O(n^2) using drsolve package.
%
% Input:
%
%   c,r -- first column c and first row r of the Toeplitz matrix T
%   b   -- RHS to solve for
%
% Output:
%
%   x   -- solution to Tx = b
%

n = length(c);

if n==1
    if r ~= c
        warning('Inconsistent Toeplitz data, column information wins');
    end
    % Special case not treated by drsolve's tsolve
    if c==0.0
        x = NaN;
    else
        x = b/c;
    end
    return;
end

% Pass through to drsolve Toeplitz solver.
x = tsolve(c,r,b,1);

% One step of iterative refinement.
for ii=1:1
    res = toepmult(c,r,x) - b;
    xx = tsolve(c,r,res,1);
    x = x - xx;
end
end