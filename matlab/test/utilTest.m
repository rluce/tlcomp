function tests = utilTest
% Test various utility functions

tests = functiontests(localfunctions);

end


function test_is_exact_toeplitz(testCase)
testCase.assertTrue(is_exact_toeplitz([]));
testCase.assertTrue(is_exact_toeplitz(0));
testCase.assertTrue(is_exact_toeplitz(1));
testCase.assertTrue(is_exact_toeplitz(1i));

testCase.assertTrue(is_exact_toeplitz([0,1; 1,0]));

testCase.assertTrue(is_exact_toeplitz(eye(5)));
testCase.assertTrue(is_exact_toeplitz(zeros(3)));
testCase.assertTrue(is_exact_toeplitz(ones(4)));
testCase.assertTrue(is_exact_toeplitz( diag((rand + 1i*rand) * ones(8))));

[c,r] = random_toeplitz(8,8);
testCase.assertTrue(is_exact_toeplitz(toeplitz(c,r)));

testCase.assertFalse(is_exact_toeplitz(rand(5)));
testCase.assertFalse(is_exact_toeplitz(eye(6) + eps * diag(rand(6))));
testCase.assertFalse(is_exact_toeplitz([0,0;0,1]));
testCase.assertFalse(is_exact_toeplitz([1 + eps, 1; 1, 1]));

end
