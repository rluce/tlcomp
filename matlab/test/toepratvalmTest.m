function tests = toepratvalmTest
% Test rational function evaluation

tests = functiontests(localfunctions);

end


function test_singleton(testCase)
c = -5;
r = -5;
p = 1;
q = 1;
[G,B] = toepratvalm(c,r,p,q);
testCase.assertEqual(stein_reconstruction(G,B), 1);

c = -5;
r = -5;
p = [1,0];
q = 1;
[G,B] = toepratvalm(c,r,p,q);
testCase.assertEqual(stein_reconstruction(G,B), -5);



c = -5;
r = -5;
p = [1,0];
q = 1;
[G,B] = toepratvalm(c,r,p,q);
testCase.assertEqual(stein_reconstruction(G,B), -5);

c = -5;
r = -5;
p = 1;
q = [1,0];
[G,B] = toepratvalm(c,r,p,q);
testCase.assertEqual(stein_reconstruction(G,B), -1/5, 'AbsTol', eps);



end

function test_zero(testCase)
c = 1;
r = 1;
p = 0;
q = 1;
[G,B] = toepratvalm(c,r,p,q);
testCase.assertEqual(stein_reconstruction(G,B), 0);



c = 0;
r = 0;
[G,B] = toepratvalm(c,r,p,q);
testCase.assertEqual(stein_reconstruction(G,B), 0);

p = [1,1,1,0];
q = [0.5, 1];
[G,B] = toepratvalm(c,r,p,q);
testCase.assertEqual(stein_reconstruction(G,B), 0);


n = 7;
c = zeros(n,1);
r = zeros(n,1);
[G,B] = toepratvalm(c,r,p,q);
% NOTE:  Somewhere in the transformation to Cauchy-like or inside drsolve's
% clsolve, a small perturbation results in a slightly nonzero solution.
% Maybe we can do better, but this is not priority.
testCase.verifyEqual(stein_reconstruction(G,B), zeros(n));
testCase.assertEqual(stein_reconstruction(G,B), zeros(n), 'AbsTol', 10*n*eps);



end

function test_identity(testCase)

n = 15;
e1 = zeros(n,1);
e1(1) = 1;

p = 1;
q = 1;
[G,B] = toepratvalm(e1, e1, p, q);
A = stein_reconstruction(G,B);
testCase.assertEqual(A, eye(n), 'AbsTol', eps);

p = [1,0];
q = 1;
[G,B] = toepratvalm(e1, e1, p, q);
A = stein_reconstruction(G,B);
testCase.assertEqual(A, eye(n), 'AbsTol', eps);

p = [1,0];
q = -2;
[G,B] = toepratvalm(e1, e1, p, q);
A = stein_reconstruction(G,B);
testCase.assertEqual(A, -eye(n)/2, 'AbsTol', eps);


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
[G,B] = toepratvalm(c, r, p, 1);
pT = stein_reconstruction(G, B);

testCase.assertEqual(pT, pT_true, 'AbsTol', 32*n*eps);

end

function test_inverse_complex(testCase)

n = 15;
[c,r,T] = lps_example3(n);

c = 1i * c;
r = 1i * r;
T = 1i * T;

p = 1;
q = [1, 0];

[G,B] = toepratvalm(c,r,p,q);
A = stein_reconstruction(G,B);
A_true = inv(T);
testCase.assertEqual(A, A_true, 'RelTol', 256*n*eps, 'AbsTol', 8*n*eps);

end

function test_inverse_real(testCase)

n = 15;
[c,r,T] = lps_example3(n);

p = 1;
q = [1, 0];

[G,B] = toepratvalm(c,r,p,q);
A = stein_reconstruction(G,B);
A_true = inv(T);
testCase.assertEqual(A, A_true, 'RelTol', 256*n*eps, 'AbsTol', 8*n*eps);

end



function test_small_complex(testCase)

c = [-1i; 1+1i];
r = [-1i, 4];

p = [ 1i, 1 + 1i];
q = [1, -3i, -1];

T = toeplitz(c,r);

pA = polyvalm(p, T);
qA = polyvalm(q, T);
rA_true = qA\pA;

[G, B] = toepratvalm(c,r,p,q);
rA = stein_reconstruction(G,B);
testCase.assertEqual(rA, rA_true, 'AbsTol', 32*eps);


end

function test_cayley(testCase)

n = 16;

% Hermitian Toeplitz matrix of spectral radius < 1
c = randn(n,1) + 1i * randn(n,1);
c(1) = real(c(1));
r = c';
T = toeplitz(c,r);
T = T / (1.1 * norm(T));
c = T(:,1);
r = T(1,:);

% Cayley transform
p = [1, -1i];
q = [1, +1i];

pA = polyvalm(p, T);
qA = polyvalm(q, T);
rA_true = qA\pA;

[G, B] = toepratvalm(c,r,p,q);
rA = stein_reconstruction(G,B);
testCase.assertEqual(rA, rA_true, 'RelTol', 32*n*eps, 'AbsTol', 8*n*eps);

end

function test_deg3_pade_real(testCase)

n = 18;

[~,~,T] = random_toeplitz(n,n);
T = real(T);

T = T / (2 * norm(T));
c = T(:,1);
r = T(1,:);

% Coeffs for (3,3) Pade approximant
[q,p] = padecoef(1,3);
pT = polyvalm(p, T);
qT = polyvalm(q, T);
rT_true = qT\pT;

[G, B] = toepratvalm(c,r,p,q);
rT = stein_reconstruction(G,B);
testCase.assertEqual(rT, rT_true, 'RelTol', n*eps, 'AbsTol', 8*n*eps);


end



function test_deg3_pade_complex(testCase)

n = 18;

[~,~,T] = random_toeplitz(n,n);

T = T / (2 * norm(T));
c = T(:,1);
r = T(1,:);

% Coeffs for (3,3) Pade approximant
[q,p] = padecoef(1,3);
pT = polyvalm(p, T);
qT = polyvalm(q, T);
rT_true = qT\pT;

[G, B] = toepratvalm(c,r,p,q);
rT = stein_reconstruction(G,B);
testCase.assertEqual(rT, rT_true, 'RelTol', 32*n*eps, 'AbsTol', 8*n*eps);


end
