function tests = tleyeTest
tests = functiontests(localfunctions);
end

function test_tleye(testCase)

E = tleye(0);
testCase.assertEqual(class(E), 'TLMat');
testCase.assertTrue(isempty(full(E)));

E = tleye(1);
testCase.assertEqual(class(E), 'TLMat');
testCase.assertEqual(full(E),1);

E = tleye(9);
testCase.assertEqual(class(E), 'TLMat');
testCase.assertEqual(full(E),eye(9));
testCase.assertEqual(drank(E), 1);
end