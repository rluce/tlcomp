function tests = toeplksquareTest

tests = functiontests(localfunctions);

end

function test_singleton(TestCase)

% 1x1 zero
G = 0;
B = 0;

[Gs, Bs] = toeplksquare(G, B);
T2 = toeplkreconstruct(Gs, Bs);
TestCase.assertEqual(T2, 0);


% 1x1 identity
G = 1;
B = 1;

runallalgs(TestCase, 1, G, B)


end

function test_identity(TestCase)
n=15;
e1 = zeros(n,1);
e1(1) = 1.0;
[G, B] = toepgen(e1,e1);
T2_true = eye(n);
runallalgs(TestCase, T2_true, G, B, 'AbsTol', n * eps)


end

function test_random_real(TestCase)
n = 12;
c = randn(n,1);
r = randn(n,1);
r(1) = c(1);
T = toeplitz(c,r);
T2_true = T*T;

[G, B] = toepgen(c,r);
runallalgs(TestCase, T2_true, G, B, 'RelTol', 1024 * n * eps)
end


function test_random_complex(TestCase)
n = 12;
c = randn(n,1) + 1i*randn(n,1);
r = randn(n,1) + 1i*randn(n,1);
r(1) = c(1);
T = toeplitz(c,r);
T2_true = T*T;
[G, B] = toepgen(c,r);
runallalgs(TestCase, T2_true, G, B, 'RelTol', 128 * n * eps)

end

function runallalgs(testCase, Ts_true, G, B, varargin)

for alg = {'full', 'fft'}
    [Gs, Bs] = toeplksquare(G, B, alg{:});
    Ts = toeplkreconstruct(Gs, Bs);
    testCase.assertEqual(Ts, Ts_true, varargin{:});
end

end

