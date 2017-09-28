function tests = displaceTest

tests = functiontests(localfunctions);

end


function test_empty(testCase)
D = displace([]);
testCase.assertTrue(isempty(D));

end

function test_onebyone(testCase)
A = 1;
D_true = 2;
D = displace(A);
testCase.assertEqual(D, D_true, 'absTol', eps);

A = -1i;
D_true = -2i;
D = displace(A);
testCase.assertEqual(D, D_true, 'absTol', eps);

end

function test_random(testCase)

c = randn(3,1) + 1i * randn(3,1);
r = randn(3,1) + 1i * randn(3,1);
r(1) = c(1);

T = toeplitz(c,r);
D_true = [
    c(3) - r(2), c(2) - r(3), c(1) + r(1);
    0          , 0          , c(2) + r(3);
    0          , 0          , c(3) + r(2);
    ];
D = displace(T);
testCase.assertEqual(D, D_true, 'absTol', eps);

end