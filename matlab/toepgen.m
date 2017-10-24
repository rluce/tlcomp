function [G,B] = toepgen(c, r)
% function [G, B] = toepgen(c, r)
%
% Return a generator of a Toeplitz matrix T.
%
% Input:
%   c  -- first column of T.
%   r  -- first row of T.
%
% The first element of c and r must be equal.
%
% The generator is not unique, one possibility is spelled out on p. 1564 of
%
%  Gohberg, I., Kailath, T., & Olshevsky, V. (1995). Fast Gaussian
%  elimination with partial pivoting for matrices with displacement
%  structure. Mathematics of Computation, 64(212), 1557?1576.

if c(1) ~= r(1)
    error('tlcomp:InconsistentInput', ...
        'first elements of c and r must be equal');
end

% Normalize to column vectors
c = c(:);
r = r(:);

n = length(r);

e1 = zeros(n,1);
en = zeros(n,1);
e1(1) = 1;
en(n) = 1;

G = [c + [0;r(end:-1:2)], e1];
B = conj([en, c(end:-1:1) - [r(2:end); 0]]);

end
