function tests = utilTest
% Test various utility functions

tests = functiontests(localfunctions);

end

function test_stein_generator1(testCase)
% Test corner case n=1

c = 1;
r = 1;
[G, B] = stein_generator(c,r);
G_true = [1, 1];
B_true = [1, 0];
testCase.assertEqual(G, G_true, 'AbsTol', eps);
testCase.assertEqual(B, B_true, 'AbsTol', eps);

% Complex variation
c = 1i;
r = 1i;
[G, B] = stein_generator(c,r);
G_true = [1i, 1];
B_true = [1, 0];
testCase.assertEqual(G, G_true, 'AbsTol', eps);
testCase.assertEqual(B, B_true, 'Abstol', eps);
end

function test_stein_generator2(testCase)
% Fail hard on inconsistent data
c = [1;2;3];
r = [4,5,6];
testCase.assertError( @() stein_generator(c,r), 'funmd:InconsistentInput');
end

function test_stein_generator3(testCase)
% Simple 3x3 example

c = [3, -2i, 0];
r = [3, 1i, 1+1i];
G_true = [
    3, 1;
    -2i, 0;
    0, 0;
    ];
B_true = [
    1, 0;
    0, -1i;
    0, 1-1i;
    ];

[G,B] = stein_generator(c,r);
testCase.assertEqual(G, G_true, 'AbsTol', eps);
testCase.assertEqual(B, B_true, 'AbsTol', eps);
end

function test_reconstruct(testCase)
T = gallery('prolate', 15, 0.51);
c = T(:,1);
r = T(1,:);

[G,B] = stein_generator(c,r);
T2 = stein_reconstruction(G,B);
testCase.assertEqual(T2, T, 'RelTol', 2*eps);
end

function test_compute_generator1(testCase)
% Random test
n = 9;
[~,~,T] = random_toeplitz(n,n);
[G, B] = compute_generator(T);

testCase.assertEqual([n,2], size(G));
testCase.assertEqual([n,2], size(B));

Z = downshift(n);
D = G * B';
D_true = T - Z*T*Z';
testCase.assertEqual(D, D_true, 'RelTol', 15*eps, 'AbsTol', 32*eps);
end

function test_compute_generator2(testCase)
% Test useful rank-1 detection
n = 13;
alpha = 3 - 2i;
T = alpha * eye(n);
[G, B] = compute_generator(T);
testCase.assertEqual([n,1], size(G));
testCase.assertEqual([n,1], size(B));

D_true = zeros(n);
D_true(1,1) = alpha;
D = G*B';
testCase.assertEqual(D, D_true, 'RelTol', 15*eps);

end

function test_compute_generator3(testCase)
% Test useful rank-5 detection
n = 11;
U = orth(randn(n,5) + 1i * randn(n,5));
V = orth(randn(n,5) + 1i * randn(n,5));

T = stein_reconstruction(U,V);
[G,B] = compute_generator(T);
testCase.assertEqual([n,5], size(G));
testCase.assertEqual([n,5], size(B));

D_true = U * V';
D = G*B';
testCase.assertEqual(D, D_true, 'RelTol', 15*eps, 'AbsTol', 10*eps);
end

function test_compute_generator4(testCase)
% Test truncation
n = 11;
U = orth(randn(n,3) + 1i * randn(n,3));
V = orth(randn(n,3) + 1i * randn(n,3));

U = U * diag([1,1,1e-8]);
T = stein_reconstruction(U,V);
[G,B] = compute_generator(T, 2);
testCase.assertEqual([n,2], size(G));
testCase.assertEqual([n,2], size(B));
Tapprox = stein_reconstruction(G,B);
testCase.assertTrue( norm(T - Tapprox) <= n * 1e-8 );
end

function test_compute_generator5(testCase)
% Requested rank too large, generators still must have correct size
n = 11;
U = orth(randn(n,3) + 1i * randn(n,3));
V = orth(randn(n,3) + 1i * randn(n,3));

T_true = stein_reconstruction(U,V);
[G,B] = compute_generator(T_true, 5);
testCase.assertEqual([n,5], size(G));
testCase.assertEqual([n,5], size(B));

T = stein_reconstruction(G,B);
testCase.assertEqual(T, T_true, 'RelTol', 15*eps, 'AbsTol', 10*eps);
end

function test_generator6(testCase)
% Fail hard on inconsistent data
[~,~,T] = random_toeplitz(8,8);
testCase.assertError( @() compute_generator(T,10), 'funmd:InconsistentInput');
end


function test_tleye(testCase)

E = tleye(0);
testCase.assertEqual(class(E), 'TLMatStein');
testCase.assertTrue(isempty(full(E)));

E = tleye(1);
testCase.assertEqual(class(E), 'TLMatStein');
testCase.assertEqual(full(E),1);

E = tleye(9);
testCase.assertEqual(class(E), 'TLMatStein');
testCase.assertEqual(full(E),eye(9));
testCase.assertEqual(drank(E), 1);


end
