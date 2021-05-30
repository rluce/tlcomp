function tests = gencompressTest

tests = functiontests(localfunctions);

end

function test_singleton(testCase)
G = 1;
B = 1;

[Gout, Bout, res_sv] = gencompress(G,B,1);
testCase.assertEqual(Gout, G);
testCase.assertEqual(Bout, B);
testCase.assertEqual(res_sv, 0.0);


[Gout, Bout, res_sv] = gencompress(G,B,6);
testCase.assertEqual(Gout, G);
testCase.assertEqual(Bout, B);
testCase.assertEqual(res_sv, 0.0);

end

function test_identity(testCase)
n = 9;
G = eye(n);
B = eye(n);

[Gred, Bred, res_sv] = gencompress(G,B,n);
testCase.assertEqual(size(Gred), [n,n]);
testCase.assertEqual(size(Bred), [n,n]);
testCase.assertEqual(Gred, G, 'AbsTol', n*eps);
testCase.assertEqual(Bred, B, 'AbsTol', n*eps);
testCase.assertLessThan(res_sv, n*eps);

[Gred, Bred, res_sv] = gencompress(G,B,n-1);
testCase.assertEqual(size(Gred), [n,n-1]);
testCase.assertEqual(size(Bred), [n,n-1]);
testCase.assertEqual(Gred, eye(n,n-1), 'AbsTol', n*eps);
testCase.assertEqual(Bred, eye(n,n-1), 'AbsTol', n*eps);
testCase.assertLessThan(res_sv, 1+n*eps);

[Gred, Bred, res_sv] = gencompress(G,B,1);
testCase.assertEqual(size(Gred), [n,1]);
testCase.assertEqual(size(Bred), [n,1]);
testCase.assertEqual(Gred, eye(9,1), 'AbsTol', n*eps);
testCase.assertEqual(Bred, eye(9,1), 'AbsTol', n*eps);
testCase.assertLessThan(res_sv, n-1+n*eps);


end


function test_zerocols(testCase)
n = 14;
G = [randn(n,3), zeros(n,3)];
B = [randn(n,3), zeros(n,3)];
G = G(:,randperm(6));
B = B(:,randperm(6));

[Gred, Bred, res_sv] = gencompress(G,B,3);
testCase.assertEqual(size(Gred), [n,3]);
testCase.assertEqual(size(Bred), [n,3]);
nfact = 2*norm(G*B');
testCase.assertEqual(Gred * Bred', G * B', ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertLessThan(res_sv, 16*eps);

end

function test_truncation_error(testCase)
n = 7;
S = diag([3,2,1,0]);
G = orth(randn(n,4));
B = orth(randn(n,4));

B = B*S;

[~, ~, res_sv] = gencompress(G,B,4);
testCase.assertLessThan(res_sv, 5*eps);

[~, ~, res_sv] = gencompress(G,B,3);
testCase.assertLessThan(res_sv, 5*eps);

[~, ~, res_sv] = gencompress(G,B,1);
testCase.assertEqual(res_sv, 2, 'RelTol', 1e-8);


end

function test_random_complex(testCase)
n = 15;
r = 3;
R = 8;

G = (randn(n,r) + 1i * randn(n,r)) * randn(r,R);
B = (randn(n,r) + 1i * randn(n,r)) * randn(r,R);

[Gred, Bred] = gencompress(G,B,r+1);
testCase.assertEqual(size(Gred), [n,r+1]);
testCase.assertEqual(size(Bred), [n,r+1]);
testCase.assertEqual(Gred*Bred', G*B', 'RelTol', 64*n*eps);

[Gred, Bred] = gencompress(G,B,r-1);
testCase.assertEqual(size(Gred), [n,r-1]);
testCase.assertEqual(size(Bred), [n,r-1]);
testCase.assertGreaterThan(norm(Gred*Bred' - G*B', 'fro'), 1e-3);

end

function test_redundant_toeplitz(testCase)

n = 13;
[c,r,T] = random_toeplitz(n,n);
[G,B] = toepgen(c,r);

% These are generators for 2 * T
GG = [G, G];
BB = [B, B];

[GGred, BBred, res_sv] = gencompress(GG, BB, 2);
testCase.assertLessThan(res_sv, 10*eps);
TT = toeplkreconstruct(GGred, BBred);
testCase.assertEqual(TT, 2*T, 'RelTol', 20*n*eps);

end

function test_rank_too_big(testCase)
G = [1,1,1,1];
B = [1,1,1,1];

% Requested length is big, we accept a shorter one.
[G, B] = gencompress(G,B,2);
testCase.assertEqual(size(G,2), 1);
testCase.assertEqual(size(B,2), 1);
T = toeplkreconstruct(G,B);
testCase.assertEqual(T, 2);


end
