function tests = vapplyTest
% Test various utility functions

tests = functiontests(localfunctions);

end

function test_arguments(testCase)
a = 1;
testCase.assertError( @() vapply(a, 'asdf'), 'expmt:InconsistentInput');
end


function test_singleton(testCase)

a = pi;

v = vapply(a);
testCase.assertEqual(v, -a, 'AbsTol', eps);

v = vapply(a, 'inv');
testCase.assertEqual(v, -a, 'AbsTol', eps);

end

function test_multiple_singleton(testCase)

a = [0.5, 1.0];
v_true = -a;

v = vapply(a);
testCase.assertEqual(v, v_true);

v = vapply(a, 'inv');
testCase.assertEqual(v, v_true);
end


function test_unitvec(testCase)

n = 349;
x = zeros(n,1);
x(1) = 1.0;

v_true = zeros(n,1);
v_true(1) = -1;
v_true(2) = 1;

v = vapply(x);
testCase.assertEqual(v, v_true, 'AbsTol', eps);

v = vapply(x, 'inv');
testCase.assertEqual(v, -ones(n,1), 'AbsTol', eps');


end


function test_onesvec(testCase)

n = 35;
x = ones(n,1);

v_true = zeros(n,1);
v_true(1) = -1.0;

v = vapply(x);
testCase.assertEqual(v, v_true, 'AbsTol', eps);

v_true = - (1:n)';
v = vapply(x, 'inv');
testCase.assertEqual(v, v_true, 'AbsTol', eps);
end


function test_10randoms(testCase)

n = 15;
V = downshift(n) - eye(n);
Vinv = -tril(ones(n));

num_tries = 10;

for k=1:num_tries
    x = randn(15,1);
    
    v_true = V*x;
    v = vapply(x);
    testCase.assertEqual(v, v_true, 'AbsTol', 2*n*eps);

    v_true = Vinv*x;
    v = vapply(x, 'inv');
    testCase.assertEqual(v, v_true, 'AbsTol', 2*n*eps);
end
end

function test_fourones(testCase)

n = 11;
x = ones(n,4);
y_true = zeros(n,4);
y_true(1,:) = -1;
y = vapply(x);
testCase.assertEqual(y, y_true, 'AbsTol', n*eps);

y_true = repmat(-(1:n)', 1, 4);
y = vapply(x, 'inv');
testCase.assertEqual(y, y_true, 'AbsTol', n*eps);

end

function test_fourrandoms(testCase)
n = 13;
x = randn(n,4);
V = downshift(n) - eye(n);
Vinv = - tril(ones(n));


y_true = V*x;
y = vapply(x);
testCase.assertEqual(y, y_true, 'RelTol', n*eps);

y_true = Vinv*x;
y = vapply(x, 'inv');
testCase.assertEqual(y, y_true, 'RelTol', 8*n*eps);
end
