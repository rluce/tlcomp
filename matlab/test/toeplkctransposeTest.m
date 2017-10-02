function tests = toeplkctransposeTest

tests = functiontests(localfunctions);

end

function test_empty(testCase)
G = [];
B = [];
[Gc, Bc] = toeplkctranspose(G, B);
testCase.assertTrue(isempty(Gc));
testCase.assertTrue(isempty(Bc));
end

function test_singleton(testCase)
G = 2;
B = 1;
Tc_true = 1;

[Gc, Bc] = toeplkctranspose(G, B);
Tc = toeplkreconstruct(Gc, Bc);
testCase.assertEqual(Tc, Tc_true, 'AbsTol', eps, 'RelTol', eps);

% Zero
G = 0;
B = 1;
Tc_true = 0;

[Gc, Bc] = toeplkctranspose(G, B);
Tc = toeplkreconstruct(Gc, Bc);
testCase.assertEqual(Tc, Tc_true, 'AbsTol', eps, 'RelTol', eps);

% Complex
G = 2i;
B = 1i;
Tc_true = -1i;

[Gc, Bc] = toeplkctranspose(G, B);
Tc = toeplkreconstruct(Gc, Bc);
testCase.assertEqual(Tc, Tc_true, 'AbsTol', eps, 'RelTol', eps);
end

function test_identity(testCase)
n = 6;
G = zeros(n,1);
B = zeros(n,1);
G(1) = 2;
B(1) = 1;
Tc_true = eye(n);

[Gc, Bc] = toeplkctranspose(G, B);
Tc = toeplkreconstruct(Gc, Bc);
testCase.assertEqual(Tc, Tc_true, 'AbsTol', eps, 'RelTol', eps);
end

function test_random_real(testCase)

n = 9;
d = 3;

G = randn(n,d);
B = randn(n,d);
Tc_true = toeplkreconstruct(G, B)';

[Gc, Bc] = toeplkctranspose(G, B);
Tc = toeplkreconstruct(Gc, Bc);
testCase.assertEqual(Tc, Tc_true, 'AbsTol', eps, 'RelTol', eps);
end



function test_random_complex(testCase)
n = 9;
d = 1;

G = randn(n,d) + 1i * randn(n,d);
B = randn(n,d) + 1i * randn(n,d);
Tc_true = toeplkreconstruct(G, B)';

[Gc, Bc] = toeplkctranspose(G, B);
Tc = toeplkreconstruct(Gc, Bc);
testCase.assertEqual(Tc, Tc_true, 'AbsTol', eps, 'RelTol', eps);
end