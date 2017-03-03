function tests = TLMatSteinTest

tests = functiontests(localfunctions);

end


function test_empty(testCase)
TL = TLMatStein([]);
testCase.assertTrue(isempty(TL.G));
testCase.assertTrue(isempty(TL.B));

TL = TLMatStein([],[]);
testCase.assertTrue(isempty(TL.G));
testCase.assertTrue(isempty(TL.B));

end