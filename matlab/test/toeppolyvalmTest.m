function tests = toeppolyvalmTest

tests = functiontests(localfunctions);

end



function test_empty_poly(testCase)

n = 14;
[c,r] = random_toeplitz(n,n);
p = [];

allalgstrial(testCase, c, r, p, zeros(n))

end

function test_zero_poly(testCase)

n = 15;

[c,r] = random_toeplitz(n,n);
p = 0;

allalgstrial(testCase, c, r, p, zeros(n))

end

function test_symmetric(testCase)

c = [1, 3/4, 1/2, 1/4, 1/5]';
r = c';
T = toeplitz(c,r);

p = [0.1, 0.5, 0.75, 1];
pT_true = polyvalm(p, T);

allalgstrial(testCase, c, r, p, pT_true, 'RelTol', 16*eps)

end



function test_singleton(testCase)

c = 0.5;
r = 0.5;
p = [1,1,1,1,1,1];

pT_true = 1.96875;

allalgstrial(testCase, c, r, p, pT_true, 'AbsTol', 8*eps)

end

function test_taylorexp(testCase)

n = 9;

% Taylor polynomial for the exponential
p = [1./5040, 1./720, 1./120, 1./24, 1./6, 1./2, 1, 1];

[c,r,T] = random_toeplitz(n,n);

% Scale down so that taylor is a good approximation -- for fun.
s = 1. / (10. * norm(T));
c = s*c;
r = s*r;
T = s*T;

pT_true = expm(T);

allalgstrial(testCase, c, r, p, pT_true, 'AbsTol', 32*n*eps)

end

function test_linear_poly_complex(testCase)

n = 41;
[c,r,T] = random_toeplitz(n,n);
p = [3.2 - 0.5i, 2*1i];

pT_true = polyvalm(p, T);

allalgstrial(testCase, c, r, p, pT_true, 'AbsTol', 32*n*eps)


end

function test_degen_poly_complex(testCase)

n = 41;
[c,r,T] = random_toeplitz(n,n);
p = [0, 3.2 - 0.5i, 2*1i];

pT_true = polyvalm(p, T);

allalgstrial(testCase, c, r, p, pT_true, 'AbsTol', 32*n*eps)

end


function test_linear_poly_real(testCase)

n = 41;
c = randn(n,1);
r = randn(n,1);
r(1) = c(1);
T = toeplitz(c,r);

p = [-1.77, 3.211];

pT_true = polyvalm(p, T);

allalgstrial(testCase, c, r, p, pT_true, 'AbsTol', 32*n*eps)

end


function test_degen_poly_real(testCase)

n = 41;
c = randn(n,1);
r = randn(n,1);
r(1) = c(1);
T = toeplitz(c,r);

p = [0, -1.77, 3.211];

pT_true = polyvalm(p, T);

allalgstrial(testCase, c, r, p, pT_true, 'AbsTol', 32*n*eps)

end

function test_degen_poly_simple_real(testCase)

n = 4;
c = [1,0,0,0];
r = [1,1,0,0];
T = toeplitz(c,r);

p = [0, 1, 1];

pT_true = polyvalm(p, T);

allalgstrial(testCase, c, r, p, pT_true, 'AbsTol', 32*n*eps)

end


function test_quadpoly_real(testCase)

c = [1,2,3,4];
r = [1,2,3,4];

p = [1,1,1];

pT_true = polyvalm(p, toeplitz(c,r));

allalgstrial(testCase, c, r, p, pT_true, 'RelTol', 100*eps)


end

function test_quadpoly_complex(testCase)

c = [1,0,0];
r = [1,1,0];

p = 1i * [1,1,1];

pT_true = polyvalm(p, toeplitz(c,r));

allalgstrial(testCase, c, r, p, pT_true, 'RelTol', 100*eps)


end

function test_tridiagonal_matrix(testCase)

c = [-1 + 2i, 1, 0, 0];
r = c;

p = 1i * [-1i, 1 - 1i, 2];

pT_true = polyvalm(p, toeplitz(c,r));

allalgstrial(testCase, c, r, p, pT_true, 'RelTol', 100*eps, 'AbsTol', 4*eps)


end

function allalgstrial(testCase, c, r, p, pT_true, varargin)

for alg = {'full', 'reduced'};%, 'horner'}
    [G, B] = toeppolyvalm(c, r, p, alg{:});
    pT = stein_reconstruction(G, B);
    testCase.assertEqual(pT, pT_true, varargin{:});
end

end
