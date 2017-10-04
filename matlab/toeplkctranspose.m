function [Gt, Bt] = toeplkctranspose(G, B)
% [Gt, Bt] = toeplkctranspose(G, B)  -- conjugate transepose
%
% Given a generator (G, B) of a TL matrix T, a generator for T' is
% returned.
%
% Input:
%
%   G, B -- generator for TL matrix T
%
% Output:
%
%   Gt, Bt -- generator for TL matrix T'
%
% NOTE: This operation increases the generator length by 2, and the
% displacement rank may also increase by 2.  Thus the obtained generator
% should be recompressed after the transposition, if minimality is
% important.
%
% See also: gencompress.m

% TODO:  We could be a bit smarter here and detect whether G has a column
% with nonzero pattern of e1, and, similar, whether B has a column of
% nonzero pattern e_n.  In that case we can merge the extraneous columns in
% the generator matrices.  How about A*e1 and A'*en?

n = size(G, 1);

if n == 0
    % Corner case, but preserve the shape of the empty matrices.
    Gt = B;
    Bt = G;
    return;
end

% First/last unit vector
e1 = zeros(n,1);
en = zeros(n,1);
e1(1) = 1;
en(n) = 1;

% For generator transformation
Zp = fcirculant(n,1);
Zm = fcirculant(n,-1);
assert(issparse(Zp) && issparse(Zm));

% First/last column of A/A'
A_e1 = toeplkmult(G, B, e1, false, 'fft');
At_en = toeplkmult(G, B, en, true, 'fft');

% Generator formulas.  Can be found from manipulating the displacement
% equation
%
%   Zp * A - A * Zm = G*B' <=> A' * Zp' - Zm' * A' = B*G' <=> ...
%
% until an expression for the Zp/Zm displacement of A' is obtained.
Gt = Zp  * [B, At_en, en];
Bt = Zm' * [G, -2*e1, -2 * A_e1];

% If n==1 then the obtained generators could be sparse.
Gt = full(Gt);
Bt = full(Bt);

end
