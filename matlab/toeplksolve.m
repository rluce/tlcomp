function sol = toeplksolve(G, B, rhs)
% sol = toeplksolve(G, B, rhs)
%
% Solve a Toeplitz-like system of equations.
%
% Input:
%   G, B   -- the Stein-Z generator of the system matrix
%   rhs    -- right hand side of the system (multiple columns allowed)
%
% Output:
%   sol    -- solution to the system

% Use partial pivoting by default

if size(G,1) == 1
    % Special code path for 1x1 matrices
    T = G*B';
    if T == 0
        % system may be inconsistent, FIXME
        sol = zeros(size(rhs));
        return;
    end
    sol = rhs / T;
    return;
end


piv = 4;
[GG, BB] = stein2sylv(G, B);

% FIXME is this really necessary?
[GG, BB] = gencompress(GG, BB);
sol = tlsolve(GG, BB, rhs, piv);
res = toeplkmult(G, B, sol) - rhs;
cor = tlsolve(GG, BB, res, piv);
sol = sol - cor;

end