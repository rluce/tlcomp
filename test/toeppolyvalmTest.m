function tests = toeppolyvalmTest

tests = functiontests(localfunctions);

end



function test_empty_poly(testCase)
% Matlab convention: empty poly of a matrix is the zero matrix of
% compatible size.
n = 14;
[c,r] = random_toeplitz(n,n);
pT_true = zeros(n);

[G, B] = toeppolyvalm(c, r, []);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true);

end

function test_zero_poly(testCase)

n = 15;

[c,r] = random_toeplitz(n,n);
p = 0;
pT_true = zeros(n);

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true);
end

function test_deg_zero(testCase)

n = 7;
c = randn(n,1);
r = randn(n,1);

E = eye(n);

for p = [0, 1, 1i, -1, -1i, rand + 1i * rand]
    [G, B] = toeppolyvalm(c, r, p);
    pT = toeplkreconstruct(G, B);
    testCase.assertEqual(pT, p * E);
end

end

function test_deg_one(testCase)

n = 8;
c = randn(n, 1) + 1i * randn(n, 1);
r = randn(n, 1) + 1i * randn(n, 1);
c(1) = r(1);
E = eye(n);
T = toeplitz(c, r);
nT = norm(T);

for p0 = [0, 1, 1i, -1, -1i, randn + 1i * randn]
    for p1 = [0, 1, 1i, -1, -1i, randn + 1i * randn]
        [G, B] = toeppolyvalm(c, r, [p1, p0]);
        pT = toeplkreconstruct(G, B);
        testCase.assertEqual(pT, p0 * E + p1 * T, ...
            'AbsTol', 4*nT*eps, 'RelTol', 4*nT*eps);
    end
end
end

function test_symmetric(testCase)

c = [1, 3/4, 1/2, 1/4, 1/5]';
r = c';
T = toeplitz(c,r);

p = [0.1, 0.5, 0.75, 1];
pT_true = polyvalm(p, T);

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', 8*eps, 'RelTol', 8*eps);
end



function test_singleton(testCase)

c = 0.5;
r = 0.5;
p = [1,1,1,1,1,1];

pT_true = 1.96875;

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true);

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

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', 512*eps, 'RelTol', 512*eps);

end

function test_linear_poly_complex(testCase)

n = 41;
[c,r,T] = random_toeplitz(n,n);
nT = norm(T);
p = [3.2 - 0.5i, 2*1i];

pT_true = polyvalm(p, T);

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', nT*4*eps, 'RelTol', nT*4*eps);

end

function test_degen_poly_complex(testCase)

n = 41;
[c,r,T] = random_toeplitz(n,n);
nT = norm(T);
p = [0, 3.2 - 0.5i, 2*1i];

pT_true = polyvalm(p, T);

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', nT*16*eps, 'RelTol', nT*16*eps);
end


function test_linear_poly_real(testCase)

n = 41;
c = randn(n,1);
r = randn(n,1);
r(1) = c(1);
T = toeplitz(c,r);
nT = norm(T);

p = [-1.77, 3.211];

pT_true = polyvalm(p, T);

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', nT*4*eps, 'RelTol', nT*4*eps);
end


function test_degen_poly_real(testCase)

n = 41;
c = randn(n,1);
r = randn(n,1);
r(1) = c(1);
T = toeplitz(c,r);
nT = norm(T);

p = [0, -1.77, 3.211];

pT_true = polyvalm(p, T);

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', nT*4*eps, 'RelTol', nT*4*eps);
end

function test_degen_poly_simple_real(testCase)

n = 4;
c = [1,0,0,0];
r = [1,1,0,0];
T = toeplitz(c,r);

p = [0, 1, 1];

pT_true = polyvalm(p, T);

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', 4*eps, 'RelTol', 4*eps);
end


function test_quadpoly_real(testCase)

c = [1,2,3,4];
r = [1,2,3,4];

p = [1,1,1];

pT_true = polyvalm(p, toeplitz(c,r));

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', 128*eps, 'RelTol', 128*eps);

end

function test_quadpoly_complex(testCase)

c = [1,0,0];
r = [1,1,0];

p = 1i * [1,1,1];

pT_true = polyvalm(p, toeplitz(c,r));

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', 4*eps, 'RelTol', 4*eps);

end

function test_tridiagonal_matrix(testCase)

c = [-1 + 2i, 1, 0, 0];
r = c;

p = 1i * [-1i, 1 - 1i, 2];

pT_true = polyvalm(p, toeplitz(c,r));

[G, B] = toeppolyvalm(c, r, p);
pT = toeplkreconstruct(G, B);
testCase.assertEqual(pT, pT_true, 'AbsTol', 16*eps, 'RelTol', 16*eps);

end