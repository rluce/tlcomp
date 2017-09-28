function tests = toeplksolvetoeplkTest
% Test toeplk \ toeplk -> toeplk

tests = functiontests(localfunctions);

end

function test_singleton(testCase)

G1 = 1;
B1 = 1;

G2 = -3;
B2 = 3;

[Gsol, Bsol] = toeplksolvetoeplk(G1, B1, G2, B2);
sol = toeplkreconstruct(Gsol, Bsol);
testCase.assertEqual(sol, -9);


end

function test_identity(testCase)

e1 = zeros(8,1);
e1(1) = 1;

G1 = e1;
B1 = e1;
G2 = e1;
B2 = e1;

[Gsol, Bsol] = toeplksolvetoeplk(G1, B1, G2, B2);
sol = toeplkreconstruct(Gsol, Bsol);
testCase.assertEqual(sol, eye(8));
end

function test_random(testCase)
n = 13;
r1 = 3;
r2 = 2;

G1 = randn(n, r1);
G1(1,1) = G1(1,1) + 1;
B1 = randn(n, r1);
G2 = randn(n, r2);
B2 = randn(n, r2);
T1 = toeplkreconstruct(G1, B1);
T2 = toeplkreconstruct(G2, B2);

[Gsol, Bsol] = toeplksolvetoeplk(G1, B1, G2, B2);
sol = toeplkreconstruct(Gsol, Bsol);
testCase.assertEqual(sol, T1\T2, 'AbsTol', 1e4*eps, 'RelTol', 1e4*eps);

end
