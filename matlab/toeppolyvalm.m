function [G, B] = toeppolyvalm(c, r, p, alg)
% [G, B] = toeppolyvalm(c, r, p)
%
% Compute Z-Stein generators of p(T) for a polynomial of a Toeplitz matrix.
%
% Arguments
%   c      -- first column of Toeplitz matrix T
%   r      -- first row of Toeplitz matrix T
%   p      -- vector of polynomial coefficients of length N+1
%   'alg'  -- evaluation method: 'full', 'reduced'
%
% Output
%   [G, B] -- Z-Stein generator pair for the polynomial
%
%        p(T) = p(1)*T^N + p(2)*T^(N-1) + ... + P(N+1)*I
%
% The Z-displacement rank of p(T) will be at most 2N, although the
% generator may be longer than that (up to N^2 + N + 1, depending on the
% actual implementation).

if nargin < 4 || isempty(alg)
    alg = 'full';
end

n = length(r);
e1 = zeros(n,1);
e1(1) = 1;

N = length(p) - 1;

if N == -1
    G = zeros(n,1);
    B = zeros(n,1);
    return;
end

if N == 0
    G = p(1) * e1;
    B = e1;
    return;
end

switch alg
    case 'full'
        [G, B] = monomial_eval(c, r, p, alg);
    case 'reduced'
        [G, B] = monomial_eval(c, r, p, alg);
    case 'horner'
        [G, B] = horner_eval(c, r, p);
    otherwise
        error('Invalid choice for alg parameter');
end

end



function [G,B] = horner_eval(c, r, p)

% Use Horner like scheme.  Computations are such that G is universal and
% the polynomial coefficients go only into B.  A priori it is not clear
% what the most stable / most efficient way of accumulation is, could be
% investigated in the future.

c = c(:);
r = r(:);
assert(c(1) == r(1));

N = length(p) - 1;
% The constant polynomial case should have been resolved by caller.
assert(N >= 1)

n = length(c);
e1 = zeros(n,1);
e1(1) = 1;
%rt = conj(r);
rt = conj(r);
rt(1) = 0;

% Generator for Horner Polynomial T1 of degree 1
G = [c, e1];
B = [conj(p(1)) * e1, [conj(p(2));conj(p(1)) * rt(2:end)]];
if (N <= 1)
    return;
end

% xi = V * T * inv(V) * e1 - c  in closed form expression
xi = [sum(r(2:end)); - r(end:-1:2)];

% Gupdt carries all powered up vectors that need to be appended
%
% NOTE: the computaiton could be arranged such that only two instead of
% three vectors need to be transformed.  However, for readibility it is
% nicer to stick as close to the writeup for this evaluation scheme.
Gupdt = [xi, e1, c];

% Initial state of the generator matrix
G = [xi, e1];
Gupdt = vapply(Gupdt, 'inv');
for k = 2:N
    Gupdt = toepmult(c,r,Gupdt);
    G = [G, vapply(Gupdt(:,1:2))]; %#ok<AGROW>
end

% Replace the second last column for the final term
G(:,end-1) = vapply(Gupdt(:,3));

% Compute B part of generator

% Bupdt carries forward products of the k-th Horner polynomial with
% (transformed) vectors e1 and tilde(r).
% Initial state: T_0' * x, where T_0 is the 0-th Horner polynomial
% p(1) * I
Bupdt = [e1, rt];
Bupdt = conj(p(1)) * vapply(Bupdt, 'inv');
for k=2:N
    % Multiplication T_k' * x = T' * (T_{k-1}' * x ) + conj(a_k) * I
    Bupdt = toepmult(conj(r), conj(c), Bupdt)  + conj(p(k)) * Bupdt;
    B = [-vapply(Bupdt(:,1)), vapply(Bupdt(:,2)) + conj(p(k+1) * e1), B]; %#ok<AGROW>
end

end % of function horner_eval

function [G,B] = monomial_eval(c,r,p,alg)

N = length(p) - 1;

[Gpow , Bpow] = toeppowers(c, r, N, alg);

for i = 1:N
    Gpow{i} = p(N-i+1) * Gpow{i};
end
G = cat(2, Gpow{1:end});
% Magic: Add Identity shift to the generator of T^1
G(1,1) = G(1,1) + p(N+1);
B = cat(2,Bpow{1:end});

% Truncate to known rank bound
[G, B] = gencompress(G,B,2*N);

end
