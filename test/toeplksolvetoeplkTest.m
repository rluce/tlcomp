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

n = 8;
e1 = zeros(n,1);
e1(1) = 1;
en = zeros(n,1);
en(n) = 1;

G1 = e1;
B1 = en;
G2 = e1;
B2 = en;

[Gsol, Bsol] = toeplksolvetoeplk(G1, B1, G2, B2);
sol = toeplkreconstruct(Gsol, Bsol);
testCase.assertEqual(sol, eye(8), 'AbsTol', 4*eps, 'RelTol', 4*eps);
end

function test_rank1(testCase)
n = 7;
G1 = randn(n,1);
B1 = randn(n,1);
G2 = randn(n,1);
B2 = randn(n,1);

T1 = toeplkreconstruct(G1, B1);
T2 = toeplkreconstruct(G2, B2);
cT1 = cond(T1);
nT2 = cond(T2);


[Gsol, Bsol] = toeplksolvetoeplk(G1, B1, G2, B2);
sol = toeplkreconstruct(Gsol, Bsol);
testCase.assertEqual(sol, T1\T2, 'AbsTol', 4*cT1*nT2*eps, 'RelTol', 4*cT1*nT2*eps);

G1 = randn(n,1) + 1i * randn(n,1);
B1 = randn(n,1) + 1i * randn(n,1);
G2 = randn(n,1) + 1i * randn(n,1);
B2 = randn(n,1) + 1i * randn(n,1);

T1 = toeplkreconstruct(G1, B1);
T2 = toeplkreconstruct(G2, B2);
cT1 = cond(T1);
nT2 = cond(T2);

[Gsol, Bsol] = toeplksolvetoeplk(G1, B1, G2, B2);
sol = toeplkreconstruct(Gsol, Bsol);
testCase.assertEqual(sol, T1\T2, 'AbsTol', 4*cT1*nT2*eps, 'RelTol', 4*cT1*nT2*eps);


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
cT1 = cond(T1);
nT2 = norm(T2);

[Gsol, Bsol] = toeplksolvetoeplk(G1, B1, G2, B2);
sol = toeplkreconstruct(Gsol, Bsol);
epsmult = 16*cT1*nT2;
testCase.assertEqual(sol, T1\T2, 'AbsTol', epsmult*eps, 'RelTol', epsmult*eps);

end
