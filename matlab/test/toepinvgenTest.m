function tests = toepinvgenTest
% Test various utility functions

tests = functiontests(localfunctions);

end

function test_singleton(testCase)

c = 1;
r = 1;

[Ginv, Binv] = toepinvgen(c,r);

Tinv = toeplkreconstruct(Ginv, Binv);
testCase.assertEqual(Tinv, 1, 'AbsTol', eps);


c = -1;
r = -1;

[Ginv, Binv] = toepinvgen(c,r);

Tinv = toeplkreconstruct(Ginv, Binv);
testCase.assertEqual(Tinv, -1, 'AbsTol', eps);


c = 1i;
r = 1i;

[Ginv, Binv] = toepinvgen(c,r);

Tinv = toeplkreconstruct(Ginv, Binv);
testCase.assertEqual(Tinv, -1i, 'AbsTol', eps);


end

function test_identity(testCase)

n = 12;
c = zeros(n,1);
r = zeros(n,1);

c(1) = 1.0;
r(1) = 1.0;

[Ginv, Binv] = toepinvgen(c,r);
I = toeplkreconstruct(Ginv, Binv);
testCase.assertEqual(I, eye(n), 'AbsTol', eps);

end


function test_random_real(testCase)

n = 7;
c = rand(n,1);
r = rand(n,1);
c(1) = r(1);

[Ginv, Binv] = toepinvgen(c,r);
testCase.assertEqual(size(Ginv), [n,2]);
testCase.assertEqual(size(Binv), [n,2]);

Tinv = toeplkreconstruct(Ginv, Binv);
T = toeplitz(c,r);
resnorm = norm(Tinv * T - eye(n), 'fro') / norm(T, 'fro');
testCase.assertLessThan(resnorm, 1e-8);

end

function test_random_complex(testCase)
n = 9;
[c,r,T] = random_toeplitz(n,n);
assert(~isreal(c(1)));

[Ginv, Binv] = toepinvgen(c,r);
testCase.assertEqual(size(Ginv), [n,2]);
testCase.assertEqual(size(Binv), [n,2]);

Tinv = toeplkreconstruct(Ginv, Binv);
resnorm = norm(Tinv * T - eye(n), 'fro') / norm(T, 'fro');
testCase.assertLessThan(resnorm, 5e-12);

end
