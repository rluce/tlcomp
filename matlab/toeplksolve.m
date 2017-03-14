function sol = toeplksolve(G, B, rhs)
% sol = toeplksolve(G, B, rhs)
%
% Solve a Toeplitz-like system of equations with dense RHS.
%
% Input:
%   G, B   -- the Stein-Z generator of the system matrix
%   rhs    -- right hand side of the system (multiple columns allowed)
%
% Output:
%   sol    -- solution to the system

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

% Use Gu pivoting w/ generator reorthogonalization
piv = 4;

% Transform to Z(1) / Z(-1) Sylvester displacement equation
[GG, BB] = stein2sylv(G, B);
% FIXME is this compression really necessary?
[GG, BB] = gencompress(GG, BB);

% Solve sytem
sol = tlsolve(GG, BB, rhs, piv);

% One step of refinement
res = toeplkmult(G, B, sol) - rhs;
cor = tlsolve(GG, BB, res, piv);
sol = sol - cor;

end