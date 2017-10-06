function tests = toeplksolveTest
% Test various utility functions

tests = functiontests(localfunctions);

end

function test_zerosolve(testCase)
G = 0;
B = 0;

b = 0;
x_true = 0;
x = toeplksolve(G, B, b);
testCase.assertEqual(x, x_true);

b = [0,0,0];
x_true = [0,0,0];
x = toeplksolve(G, B, b);
testCase.assertEqual(x, x_true);

G = [0,0,0,0,0];
B = [0,0,0,0,0];
b = [0,0,0];
x_true = [0,0,0];
x = toeplksolve(G, B, b);
testCase.assertEqual(x, x_true);
end

function test_singleton(testCase)

G = 2;
B = 1;

b = 1;
x_true = 1;
x = toeplksolve(G, B, b);
testCase.assertEqual(x, x_true);

b = [1,2,3];
x_true = [1,2,3];
x = toeplksolve(G, B, b);
testCase.assertEqual(x, x_true);

G = 2i;
B = 1;

b = [1,2,3];
x_true = -1i * [1,2,3];
x = toeplksolve(G, B, b);
testCase.assertEqual(x, x_true);
end

function test_prolate_matrix(testCase)

n = 17;
T = gallery('prolate', n, 0.35);
[G, B] = toepgen(T(:,1), T(1,:));
b = [ones(n,1), linspace(0,1,n)'];
x_true = T\b;
x = toeplksolve(G,B,b);
% Condition number around 10^5
testCase.assertEqual(x, x_true, 'RelTol', 10^6 * n*eps);

end

function test_random_toeppoly(testCase)
n = 13;

[c,r,T] = random_toeplitz(n,n);

polycoef = [0.01, 0.1i, 1];
pT = polyvalm(polycoef, T);
[G, B] = toeppolyvalm(c,r,polycoef);

b = [ones(n,1), linspace(0,1,n)'];
x_true = pT\b;
x = toeplksolve(G,B,b);
testCase.assertEqual(x, x_true, 'RelTol', 32*n*eps);
end
