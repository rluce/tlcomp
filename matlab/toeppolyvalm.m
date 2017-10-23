function [G, B] = toeppolyvalm(c, r, p)
% [G, B] = toeppolyvalm(c, r, p)
%
% Compute Zp1/Zm1 generator of p(T) for a polynomial p of a Toeplitz matrix T.
%
% Arguments
%   c      -- first column of Toeplitz matrix T
%   r      -- first row of Toeplitz matrix T
%   p      -- vector of polynomial coefficients of length N+1
%
% Output
%   [G, B] -- generator pair for the polynomial
%
%        p(T) = p(1)*T^N + p(2)*T^(N-1) + ... + P(N+1)*I
%
% The displacement rank of p(T) will be at most 2N.

n = length(r);
e1 = zeros(n,1);
e1(1) = 1;

en = zeros(n,1);
en(n) = 1;

N = length(p) - 1;

if N == -1
    G = zeros(n,1);
    B = zeros(n,1);
    return;
end

if N == 0
    G = 2 * p(1) * e1;
    B = en;
    return;
end

[G, B] = monomial_eval(c, r, p);

end


function [G,B] = monomial_eval(c,r,p)

N = length(p) - 1;

[Gpow , Bpow] = toeppowers(c, r, N);

for i = 1:N
    Gpow{i} = p(N-i+1) * Gpow{i};
end
G = cat(2, Gpow{1:end});
% Magic: Add Identity shift to the generator of T^1
G(1,1) = G(1,1) + 2*p(N+1);
B = cat(2,Bpow{1:end});

% Truncate to known rank bound
[G, B] = gencompress(G,B,2*N);

end
