function tests = toeplkprodTest

tests = functiontests(localfunctions);

end

function test_eye(testCase)
G1 = 1;
B1 = 1;
G2 = 1;
B2 = 1;

[GP, BP] = toeplkprod(G1, B1, G2, B2);
testCase.assertEqual(stein_reconstruction(GP, BP), 1);

e1 = zeros(8);
e1(1) = 1;
G1 = e1;
B1 = e1;
G2 = e1;
B2 = e1;

[GP, BP] = toeplkprod(G1, B1, G2, B2);
testCase.assertEqual(stein_reconstruction(GP, BP), eye(8));
end

function test_square(testCase)

G1 = randn(12,4);
B1 = randn(12,4);
G2 = G1;
B2 = B1;
P_true = stein_reconstruction(G1, B1)^2;
[GP, BP] = toeplkprod(G1, B1, G2, B2);
P = stein_reconstruction(GP, BP);
testCase.assertEqual(P, P_true, 'AbsTol', 1000*eps, 'RelTol', 1000*eps);
end

function test_random(testCase)
G1 = randn(12,4);
B1 = randn(12,4);
G2 = 1i * randn(12,3);
B2 = 1i * randn(12,3);
P_true = stein_reconstruction(G1, B1) * stein_reconstruction(G2, B2);
[GP, BP] = toeplkprod(G1, B1, G2, B2);
P = stein_reconstruction(GP, BP);
testCase.assertEqual(P, P_true, 'AbsTol', 1000*eps, 'RelTol', 1000*eps);
end

function test_toeplitz(testCase)

[c1, r1, T1] = random_toeplitz(8,8);
[c2, r2, T2] = random_toeplitz(8,8);
[G1, B1] = stein_generator(c1,r1);
[G2, B2] = stein_generator(c2,r2);
[GP, BP] = toeplkprod(G1, B1, G2, B2);
P = stein_reconstruction(GP, BP);
testCase.assertEqual(P, T1 * T2, 'AbsTol', 100*eps, 'RelTol', 100*eps);


end
