function tests = toepeyeTest

tests = functiontests(localfunctions);

end

function test_teye(testCase)

E = toepeye(0);
testCase.assertEqual(class(E), 'ToepMat');
testCase.assertTrue(isempty(full(E)));

E = toepeye([0,0]);
testCase.assertEqual(class(E), 'ToepMat');
testCase.assertTrue(isempty(full(E)));

E = toepeye(0,0);
testCase.assertEqual(class(E), 'ToepMat');
testCase.assertTrue(isempty(full(E)));

E = toepeye(1);
testCase.assertEqual(class(E), 'ToepMat');
testCase.assertEqual(full(E),1);

E = toepeye(1,1);
testCase.assertEqual(class(E), 'ToepMat');
testCase.assertEqual(full(E),1);

E = toepeye([1,1]);
testCase.assertEqual(class(E), 'ToepMat');
testCase.assertEqual(full(E),1);

E = toepeye(9);
testCase.assertEqual(class(E), 'ToepMat');
testCase.assertEqual(full(E),eye(9));

E = toepeye(9,9);
testCase.assertEqual(class(E), 'ToepMat');
testCase.assertEqual(full(E),eye(9));

E = toepeye([9,9]);
testCase.assertEqual(class(E), 'ToepMat');
testCase.assertEqual(full(E),eye(9));

end