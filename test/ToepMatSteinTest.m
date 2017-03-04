function tests = ToepMatSteinTest

tests = functiontests(localfunctions);

end


function test_constructor(testCase)

T = ToepMatStein([], []);
testCase.assertTrue(isempty(T.c));
testCase.assertTrue(isempty(T.r));

T = ToepMatStein(-1, -1);
testCase.assertEqual(T.c, -1);
testCase.assertEqual(T.r, -1);

[c,r] = random_toeplitz(7,7);

T = ToepMatStein(c,r);
testCase.assertEqual(T.c, c);
testCase.assertEqual(T.r, r);

end