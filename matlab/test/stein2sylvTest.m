function tests = stein2sylvTest

tests = functiontests(localfunctions);

end

function D = sylvdisplacement(T)
n = size(T, 1);

Z1 = fcirculant(n,1);
Zm1 = fcirculant(n,-1);

D = Z1 * T - T * Zm1;
D = full(D);
end

function test_singleton(testCase)

Gstein = 0;
Bstein = 0;

[Gsylv, Bsylv] = stein2sylv(Gstein, Bstein);
D = Gsylv * Bsylv';
D_true = sylvdisplacement(stein_reconstruction(Gstein, Bstein));
testCase.assertEqual(D, D_true, 'RelTol', 8*eps);


Gstein = 1;
Bstein = 1;

[Gsylv, Bsylv] = stein2sylv(Gstein, Bstein);
D = Gsylv * Bsylv';
D_true = sylvdisplacement(stein_reconstruction(Gstein, Bstein));
testCase.assertEqual(D, D_true, 'RelTol', 8*eps);

Gstein = [1,1,1];
Bstein = [1,1,1];

[Gsylv, Bsylv] = stein2sylv(Gstein, Bstein);
D = Gsylv * Bsylv';
D_true = sylvdisplacement(stein_reconstruction(Gstein, Bstein));
testCase.assertEqual(D, D_true, 'RelTol', 8*eps);
end

function test_identity(testCase)
n = 13;
e1 = zeros(n,1);
e1(1) = 1;

T = eye(n);
[Gstein, Bstein] = stein_generator(e1,e1);
[Gsylv, Bsylv] = stein2sylv(Gstein, Bstein);
D = Gsylv * Bsylv';
D_true = sylvdisplacement(T);
testCase.assertEqual(D, D_true, 'RelTol', n*eps);

end

function test_zeros(testCase)

n = 17;
x = zeros(n,1);

T = zeros(n);
[Gstein, Bstein] = stein_generator(x,x);
[Gsylv, Bsylv] = stein2sylv(Gstein, Bstein);
D = Gsylv * Bsylv';
D_true = sylvdisplacement(T);
testCase.assertEqual(D, D_true, 'AbsTol', n*eps);

end

function test_random_toepltiz(testCase)

n = 12;
[c, r, T] = random_toeplitz(n, n);
[Gstein, Bstein] = stein_generator(c,r);
[Gsylv, Bsylv] = stein2sylv(Gstein, Bstein);
D = Gsylv * Bsylv';
D_true = sylvdisplacement(T);
testCase.assertEqual(D, D_true, 'AbsTol', n*eps);

end

function test_random_data(testCase)

n = 12;
r = 5;

Gstein = orth(randn(n,r));
Bstein = orth(randn(n,r));

T = stein_reconstruction(Gstein,Bstein);
[Gsylv, Bsylv] = stein2sylv(Gstein, Bstein);
D = Gsylv * Bsylv';
D_true = sylvdisplacement(T);
testCase.assertEqual(D, D_true, 'AbsTol', n*eps);

end