function tests = toepsolveTest
% Test various utility functions

tests = functiontests(localfunctions);

end


function test_singleton(testCase)

c = 1;
r = 1;
b = 1;

x = toepsolve(c,r,b);

testCase.assertEqual(x, 1.0);

b = 0;
x = toepsolve(c,r,b);

testCase.assertEqual(x, 0.0);

end

function test_toepdata(testCase)

n = 13;

c = randn(n,1);
r = randn(n,1);
r(1) = c(1);

x_true = zeros(n,1);
x_true(1) = 1.0;
b = c;

x = toepsolve(c,r,b);

testCase.assertEqual(x, x_true, 'AbsTol', 32*n*eps);

% Solve for the transpose
b = r;
x = toepsolve(r,c,b);
testCase.assertEqual(x, x_true, 'AbsTol', 16*n*eps);

end

function test_identity(testCase)

n = 51;
c = zeros(n,1);
r = zeros(n,1);
c(1) = 1;
r(1) = 1;

b = randn(n,1);
x_true = b;
x = toepsolve(c,r,b);
testCase.assertEqual(x, x_true, 'AbsTol', n*eps);


end

function test_zeromatrix(testCase)

% TODO unclear what the desired behaviour actually is, set up only for now

n = 125;
r = zeros(n,1);
c = zeros(n,1);
b = randn(n,1);

testCase.assertWarning( @() toepsolve(c,r,b), 'drsolve:clsolve:singularMatrix');

b = randn(n,4);
testCase.assertWarning( @() toepsolve(c,r,b), 'drsolve:clsolve:singularMatrix');
end

function test_multiplerhs(testCase)
n = 13;

c = randn(n,1);
r = randn(n,1);
r(1) = c(1);

x_true = zeros(n,2);
x_true(1,1) = 1.0;
x_true(n,2) = 1.0;

b = [c, r(end:-1:1)];

x = toepsolve(c,r,b);

testCase.assertEqual(x, x_true, 'AbsTol', 64*n*eps);

% Solve for the transpose
b = [r, c(end:-1:1)];
x = toepsolve(r,c,b);
testCase.assertEqual(x, x_true, 'AbsTol', 32*n*eps);
end