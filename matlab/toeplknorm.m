function val = toeplknorm(G, B, p)
% val = toepnorm(c, r, p)  --  Compute matrix norm |T|_p of a Toeplitz-like
% matrix T.
%
% Input:
%   G, B -- generator of T
%   p    -- which norm to compute, possible choices:
%             1    : Matrix-1 norm (maximum column abs sum)
%             inf  : Matrix-inf norm (maximum row abs sum)
%             'fro': Frobenius norm
%
% Output:
%   val   -- |T|_p
%
% All of these norms are computed in O(d n^2) ops and O(n) memory.
%
% The 2-norm cannot be computed as cheaply as the norms above, which is
% why it is not supported.  For this norm consider using |normest|.

% TODO: Is there a way to compute these norms in O(d * n) ops?

switch p
    case 1
        val = toep_norm_1(G, B);
    case inf
        val = toep_norm_inf(G, B);
    case 'fro'
        val = toep_norm_fro(G, B);
    otherwise
        error('tlzstein:Unsupported', ...
            'Only matrix norms 1, inf and ''fro'' are supported');
end
end

function tnorm = toep_norm_fro(G,B)
tnorm = 0.0;

    function visitor_fro(col)
        tnorm = tnorm + sum(abs(col(:)).^2);
    end

column_visitor(G, B, @visitor_fro);
tnorm = sqrt(tnorm);
end

function tnorm = toep_norm_1(G,B)
tnorm = 0.0;

    function visitor_1(col)
        this_val = max(sum(abs(col)));
        tnorm = max([this_val, tnorm]);
    end
column_visitor(G, B, @visitor_1);
end

function tnorm = toep_norm_inf(G,B)
rowabssum = zeros(size(G,1),1);

    function visitor_inf(col)
        rowabssum = rowabssum + sum(abs(col), 2);
    end

column_visitor(G, B, @visitor_inf);
tnorm = max([0; rowabssum]); % leading zero covers corner case of empty matrix
end

function column_visitor(G, B, f)
% Traverses the the columns t_j of T(G,B), and calls f(t_j) for each such
% column.
%
% TODO: This function could possibly factored out and used for other
% purposes as well.
%
% Adapted from toeplkreconstruct

if isempty(G)
    % Then T has no columns, and so all columns are processed at this
    % point.
    return;
end

n = size(G, 1);

% Transpose for more comfi notation below
B = B';

% First column is diagonal-wrap-around cumsum of G*B
s = G * B(:,1);
for k = 2:n
    t = G * B(:,k);
    tmpval = s(end) + t(1);
    s(2:end) = s(1:end-1) + t(2:end);
    s(1) = tmpval;
end

% This is the first column of T
newvals = .5 * s;

% Remaining columns arise from first by
% cumulative-wrap-around-diagonal-sums
for k = 1:n
    t = G*B(:,k);
    % tmpvec will hold the updated values for column k+1
    tmpvec(1,1) = -t(1) + newvals(end); % Wrap around value
    tmpvec(2:n,1) = -t(2:end) + newvals(1:end-1); % No wrap in rest of col
    
    % Visit column k
    f(newvals);
    
    % Iteration upkeep: Carry over values for column k+1
    newvals = tmpvec;
end

end