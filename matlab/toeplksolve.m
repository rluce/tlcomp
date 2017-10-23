function sol = toeplksolve(G, B, rhs, ctrans)
% sol = toeplksolve(G, B, rhs)
%
% Solve a Toeplitz-like system Ax = rhs of equations with dense RHS.
%
% Input:
%   G, B   -- the generator of the system matrix A
%   rhs    -- right hand side rhs of the system (multiple columns allowed)
%   ctrans -- (optional) true/false switch.  If true, solve A'x = rhs
%             instead. Default: false.
%
% Output:
%   sol    -- solution to the system

if size(G,1) == 1
    % Special code path for 1x1 matrices
    T = 0.5 * G*B'; % 1x1 recostruction formula
    if T == 0
        % FIXME: System may be inconsistent, we should issue a warning or
        % so.
        sol = zeros(size(rhs));
        return;
    end
    sol = rhs / T;
    return;
end

% Use Gu pivoting w/ generator reorthogonalization
piv = 4;

% Solve sytem
sol = tlsolve(G, B, rhs, piv);

% One step of refinement
res = toeplkmult(G, B, sol) - rhs;
cor = tlsolve(G, B, res, piv);
sol = sol - cor;
end