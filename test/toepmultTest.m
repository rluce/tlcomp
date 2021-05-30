function tests = toepmultTest

tests = functiontests(localfunctions);

end


function test_singleton(testCase)
c = 0;
r = 0;

x = 1;
y_true = 0;
y = toepmult(c,r,x);
testCase.assertEqual(y, y_true, 'AbsTol', eps);

c = 1;
r = 1;

x = 1;
y_true = 1;
y = toepmult(c,r,x);
testCase.assertEqual(y, y_true, 'AbsTol', eps);

c = 1i;
r = 1i;

x = 1i;
y_true = -1;
y = toepmult(c,r,x);
testCase.assertEqual(y, y_true, 'AbsTol', eps);
end


function test_identity(testCase)

n = 13;

e1 = zeros(n,1);
e1(1) = 1.0;

c = e1;
r = e1;

x = e1;
y_true = e1;
y = toepmult(c,r,x);
testCase.assertEqual(y, y_true, 'AbsTol', n*eps);

x = randn(n,1);
y_true = x;
y = toepmult(c,r,x);
testCase.assertEqual(y, y_true, 'AbsTol', n*eps);

x = randn(n,1) + 1i * randn(n,1);
y_true = x;
y = toepmult(c,r,x);
testCase.assertEqual(y, y_true, 'AbsTol', n*eps);


end

function test_tril_ones(testCase)
n = 15;

c = ones(n,1);
r = zeros(n,1);
r(1) = 1;
x = ones(n,1);
y_true = (1:n)';
y = toepmult(c,r,x);
testCase.assertEqual(y, y_true, 'RelTol', n*eps);


end

function test_triltoep(testCase)
n = 14;
c = randn(n,1);
r = zeros(n,1);
r(1) = c(1);
T = toeplitz(c,r);
x = ones(n,1);
y_true = T*x;
y = toepmult(c,r,x);
testCase.assertEqual(y, y_true, 'RelTol', 128*n*eps);

end

function test_multiplerhs(testCase)
n = 5;
c = randn(n,1);
r = randn(n,1);
r(1) = c(1);
x = ones(n,2);

y = toepmult(c,r,x);

testCase.assertEqual(y(:,1), y(:,2), 'AbsTol', n*eps);

end