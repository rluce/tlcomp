function [G,B] = stein_generator(c, r)
% function [G,B] = stein_generator(c, r)
%
% Return the canonical generator of a Toeplitz matrix w.r.t Stein
% displacement.

if c(1) ~= r(1)
    error('expmt:InconsistentInput', ...
        'first elements of c and r must be the same');
end

% Normalize to column vectors
c = c(:);
r = r(:);

n = length(c);
e1 = zeros(n,1);
e1(1) = 1;

G = [c, e1];
B = [e1, [0; conj(r(2:end))]];
end