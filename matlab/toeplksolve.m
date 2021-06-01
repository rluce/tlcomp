function sol = toeplksolve(G, B, rhs, ctrans)
% sol = toeplksolve(G, B, rhs)
% sol = toeplksolve(G, B, rhs, ctrans)
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

if nargin < 4 || isempty(ctrans)
    ctrans = false;
end

n = size(G,1);
if n == 1
    % Special code path for 1x1 matrices
    T = 0.5 * G*B'; % 1x1 recostruction formula
    if ctrans
        T = T';
    end
    if T == 0
        % FIXME: System may be inconsistent, we should issue a warning or
        % so.
        sol = zeros(size(rhs));
        return;
    end
    sol = rhs / T;
    return;
end

if ctrans
    % We convert the representation of A w.r.t Zp1 / Zm1 displacement to a
    % representation of A' w.r.t. Zm1 / Zp1.
    Zp1 = fcirculant(n, 1);
    Zm1 = fcirculant(n, -1);
    GG = Zm1 * B;
    BB = Zp1' * G;
    xi = -1;
    eta = 1;
else
    % We use the natural representation.
    GG = G;
    BB = B;
    xi = 1;
    eta = -1;
end

% Use Gu pivoting w/ generator reorthogonalization
piv = 4;

% Solve sytem
sol = tlsolve(GG, BB, rhs, piv, xi, eta);

% One step of refinement
do_refine = true;
if do_refine
    res = toeplkmult(G, B, sol, ctrans) - rhs;
    cor = tlsolve(GG, BB, res, piv, xi, eta);
    sol = sol - cor;
end

end
