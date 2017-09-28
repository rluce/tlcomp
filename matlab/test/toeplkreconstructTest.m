function tests = toeplkreconstructTest
% Test various utility functions

tests = functiontests(localfunctions);

end


function test_simple(testCase)
% Simple 3x3 example

c_true = [3, -2i, 0];
r_true = [3, 1i, 1+1i];
T_true = toeplitz(c_true,r_true);

% This is a Z1/Zm1 generator for T.
G = [
    3,   1;
 1-1i,   0;
    1i,   0;
    ];
B = [
    0, 1i;
    0, -1+3i;
    1, 3;
    ];

T = toeplkreconstruct(G, B);
testCase.assertEqual(T, T_true, 'AbsTol', eps);
end

function test_reconstruct(testCase)
% toepgen / toeplkreconstruct path

T = gallery('prolate', 15, 0.51);
c = T(:,1);
r = T(1,:);

[G,B] = toepgen(c,r);
T2 = toeplkreconstruct(G,B);
testCase.assertEqual(T2, T, 'RelTol', 2*eps);
end