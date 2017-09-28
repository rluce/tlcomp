function tests = toeplkgenTest
tests = functiontests(localfunctions);
end

function test_compute_generator1(testCase)
% Random test
n = 9;
[~,~,T] = random_toeplitz(n,n);
[G, B] = toeplkgen(T);

testCase.assertEqual([n,2], size(G));
testCase.assertEqual([n,2], size(B));

D = G * B';
D_true = displace(T);
testCase.assertEqual(D, D_true, 'RelTol', 15*eps, 'AbsTol', 32*eps);
end

function test_compute_generator2(testCase)
% Recognize an exact rank 1 matrix
n = 13;
alpha = 3 - 2i;
T = alpha * eye(n);
[G, B] = toeplkgen(T);
testCase.assertEqual([n,1], size(G));
testCase.assertEqual([n,1], size(B));

D_true = zeros(n);
D_true(1,n) = 2 * alpha;
D = G*B';
testCase.assertEqual(D, D_true, 'RelTol', eps);

end

function test_compute_generator3(testCase)
% rank-5
n = 11;
U = orth(randn(n,5) + 1i * randn(n,5));
V = orth(randn(n,5) + 1i * randn(n,5));

T = toeplkreconstruct(U,V);
[G,B] = toeplkgen(T);
testCase.assertEqual([n,5], size(G));
testCase.assertEqual([n,5], size(B));

D_true = U * V';
D = G*B';
testCase.assertEqual(D, D_true, 'RelTol', 15*eps, 'AbsTol', 10*eps);
end

function test_compute_generator4(testCase)
% Test pseudo-optimality of truncation
n = 11;
U = orth(randn(n,3) + 1i * randn(n,3));
V = orth(randn(n,3) + 1i * randn(n,3));

U = U * diag([1,1,1e-8]);
T = toeplkreconstruct(U,V);
[G,B] = toeplkgen(T, 2);
testCase.assertEqual([n,2], size(G));
testCase.assertEqual([n,2], size(B));
Tapprox = toeplkreconstruct(G,B);
testCase.assertTrue( norm(T - Tapprox) <= n * 1e-8 );
end

function test_compute_generator5(testCase)
% Requested rank too large, generators still must have correct size
n = 11;
U = orth(randn(n,3) + 1i * randn(n,3));
V = orth(randn(n,3) + 1i * randn(n,3));

T_true = toeplkreconstruct(U,V);
[G,B] = toeplkgen(T_true, 5);
testCase.assertEqual([n,5], size(G));
testCase.assertEqual([n,5], size(B));

T = toeplkreconstruct(G,B);
testCase.assertEqual(T, T_true, 'RelTol', 15*eps, 'AbsTol', 10*eps);
end

function test_generator6(testCase)
% Fail hard on inconsistent data
[~,~,T] = random_toeplitz(8,8);
testCase.assertError( @() toeplkgen(T,10), ...
    'tlzstein:InconsistentInput');
end
