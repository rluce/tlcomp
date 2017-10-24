function tests = toepgenTest
% Test various utility functions

tests = functiontests(localfunctions);

end

function test_onebyone(testCase)
% Test corner case n=1

c = 1;
r = 1;
[G, B] = toepgen(c,r);
G_true = [1, 1];
B_true = [1, 1];
testCase.assertEqual(G, G_true, 'AbsTol', eps);
testCase.assertEqual(B, B_true, 'AbsTol', eps);

% Complex variation
c = 1i;
r = 1i;
[G, B] = toepgen(c,r);
G_true = [1i, 1];
B_true = [1, -1i];
testCase.assertEqual(G, G_true, 'AbsTol', eps);
testCase.assertEqual(B, B_true, 'Abstol', eps);
end

function test_inconsistent(testCase)
% Fail hard on inconsistent data
c = [1;2;3];
r = [4,5,6];
testCase.assertError( @() toepgen(c,r), 'tlcomp:InconsistentInput');
end

function test_simple(testCase)
% Simple 3x3 example

c = [3, -2i, 0];
r = [3, 1i, 1+1i];
G_true = [
    3,   1;
 1-1i,   0;
    1i,   0;
    ];
B_true = [
    0, 1i;
    0, -1+3i;
    1, 3;
    ];

[G,B] = toepgen(c,r);
testCase.assertEqual(G, G_true, 'AbsTol', eps);
testCase.assertEqual(B, B_true, 'AbsTol', eps);
end
