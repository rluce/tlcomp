function tests = istoeplitzTest

tests = functiontests(localfunctions);

end


function test_istoeplitz(testCase)
testCase.assertTrue(istoeplitz([]));
testCase.assertTrue(istoeplitz(0));
testCase.assertTrue(istoeplitz(1));
testCase.assertTrue(istoeplitz(1i));

testCase.assertTrue(istoeplitz([0,1; 1,0]));

testCase.assertTrue(istoeplitz(eye(5)));
testCase.assertTrue(istoeplitz(zeros(3)));
testCase.assertTrue(istoeplitz(ones(4)));
testCase.assertTrue(istoeplitz( diag((rand + 1i*rand) * ones(8))));

[c,r] = random_toeplitz(8,8);
testCase.assertTrue(istoeplitz(toeplitz(c,r)));

testCase.assertFalse(istoeplitz(rand(5)));
testCase.assertFalse(istoeplitz(eye(6) + eps * diag(rand(6))));
testCase.assertFalse(istoeplitz([0,0;0,1]));
testCase.assertFalse(istoeplitz([1 + eps, 1; 1, 1]));

end
