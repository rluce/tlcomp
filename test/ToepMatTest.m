function tests = ToepMatTest

tests = functiontests(localfunctions);

end


function test_constructor(testCase)

T = ToepMat([], []);
testCase.assertTrue(isempty(T.c));
testCase.assertTrue(isempty(T.r));

T = ToepMat(-1, -1);
testCase.assertEqual(T.c, -1);
testCase.assertEqual(T.r, -1);

[c,r] = random_toeplitz(7,7);

T = ToepMat(c,r);
testCase.assertEqual(T.c, c(:));
testCase.assertEqual(T.r, r(:));

c = 6;
r = 1;
testCase.assertError( @() ToepMat(c,r), 'tlzstein:InconsistentInput');

c = rand(9,1);
r = c;
r(1) = r(1) + 1e-8;
testCase.assertError( @() ToepMat(c,r), 'tlzstein:InconsistentInput');

% For the moment, we only support square matrices
c = ones(4,1);
r = ones(3,1);
testCase.assertError( @() ToepMat(c,r), 'tlzstein:InconsistentInput');

% Only vectors allowed
c = ones(4,2);
r = ones(4,2);
testCase.assertError( @() ToepMat(c,r), 'tlzstein:InconsistentInput');

end

function test_size(testCase)

T = ToepMat([], []);
s = size(T);
testCase.assertEqual(s, [0,0]);
[s1, s2] = size(T);
testCase.assertEqual(s1, 0);
testCase.assertEqual(s2, 0);

T = ToepMat(1,1);
s = size(T);
testCase.assertEqual(s, [1,1]);
[s1, s2] = size(T);
testCase.assertEqual(s1, 1);
testCase.assertEqual(s2, 1);

T = ToepMat([1,2,3], [1,0, -1]);
s = size(T);
testCase.assertEqual(s, [3,3]);
[s1, s2] = size(T);
testCase.assertEqual(s1, 3);
testCase.assertEqual(s2, 3);


end

function test_full(testCase)
T = ToepMat([],[]);
testCase.assertEqual(full(T), []);

T = ToepMat(-1,-1);
testCase.assertEqual(full(T), -1);

T = ToepMat([1,0,0], [1,0,0]);
testCase.assertEqual(full(T), eye(3));

[c,r,T] = random_toeplitz(8,8);
TM = ToepMat(c,r);
testCase.assertEqual(full(TM), T);
end

