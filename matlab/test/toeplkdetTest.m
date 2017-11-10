function tests = toeplkdetTest
tests = functiontests(localfunctions);
end

function test_singleton(testCase)

d = toeplkdet(0,0);
testCase.assertEqual(d, 0.0);

% 1x1 eye matrix
d = toeplkdet(1,2);
testCase.assertEqual(d, 1.0);

d = toeplkmult(1,pi);
testCase.assertEqual(d, pi/2);
end

function test_zero(testCase)
d = toeplkdet(zeros(8,3), zeros(8,3));
testCase.assertEqual(d, 0.0);
end

function test_identity(testCase)
[G, B] = toeplkgen(eye(9));
d = toeplkdet(G, B);
testCase.assertEqual(d, 1.0);
end

function test_random_real(testCase)
A = randn(5);
[G, B] = toeplkgen(A);
d = toeplkdet(G,B);
d_true = det(A);
testCase.assertEqual(d, d_true);
end

function test_random_complex(testCase)
A = randn(5) + 1i * rand(5);
[G, B] = toeplkgen(A);
d = toeplkdet(G,B);
d_true = det(A);
testCase.assertEqual(d, d_true);
end
