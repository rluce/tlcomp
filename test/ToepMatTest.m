function tests = ToepMatTest

tests = functiontests(localfunctions);

end


function test_constructor(testCase)

T = ToepMat([], []);
testCase.assertTrue(isempty(T.c));
testCase.assertTrue(isempty(T.r));

T = ToepMat(-1, -1);
testCase.assertEqual(T.c, -1);
testCase.assertEqual(T.r, -1);

[c,r] = random_toeplitz(7,7);

T = ToepMat(c,r);
testCase.assertEqual(T.c, c(:));
testCase.assertEqual(T.r, r(:));

c = 6;
r = 1;
testCase.assertError( @() ToepMat(c,r), 'tlzstein:InconsistentInput');

c = rand(9,1);
r = c;
r(1) = r(1) + 1e-8;
testCase.assertError( @() ToepMat(c,r), 'tlzstein:InconsistentInput');

% For the moment, we only support square matrices
c = ones(4,1);
r = ones(3,1);
testCase.assertError( @() ToepMat(c,r), 'tlzstein:InconsistentInput');

% Only vectors allowed
c = ones(4,2);
r = ones(4,2);
testCase.assertError( @() ToepMat(c,r), 'tlzstein:InconsistentInput');

end

function test_size(testCase)

T = ToepMat([], []);
s = size(T);
testCase.assertEqual(s, [0,0]);
[s1, s2] = size(T);
testCase.assertEqual(s1, 0);
testCase.assertEqual(s2, 0);

T = ToepMat(1,1);
s = size(T);
testCase.assertEqual(s, [1,1]);
[s1, s2] = size(T);
testCase.assertEqual(s1, 1);
testCase.assertEqual(s2, 1);

T = ToepMat([1,2,3], [1,0, -1]);
s = size(T);
testCase.assertEqual(s, [3,3]);
[s1, s2] = size(T);
testCase.assertEqual(s1, 3);
testCase.assertEqual(s2, 3);


end

function test_full(testCase)
T = ToepMat([],[]);
testCase.assertEqual(full(T), []);

T = ToepMat(-1,-1);
testCase.assertEqual(full(T), -1);

T = ToepMat([1,0,0], [1,0,0]);
testCase.assertEqual(full(T), eye(3));

[c,r,T] = random_toeplitz(8,8);
TM = ToepMat(c,r);
testCase.assertEqual(full(TM), T);
end

function test_add_scalar(testCase)

T = ToepMat([],[]);
testCase.assertEqual(class(T + exp(1)), 'ToepMat');
testCase.assertTrue(isempty(full(T+exp(1))));
testCase.assertEqual(class(exp(1) + T), 'ToepMat');
testCase.assertTrue(isempty(full(exp(1) + T)));
testCase.assertEqual(class(exp(1) - T), 'ToepMat');
testCase.assertTrue(isempty(full(exp(1) - T)));

T = ToepMat(2,2);
T = T - 1;
testCase.assertEqual(T.c, 1);
testCase.assertEqual(T.r, 1);

T = T + 0;
testCase.assertEqual(T.c, 1);
testCase.assertEqual(T.r, 1);

[c, r, T] = random_toeplitz(5,5);
TM = ToepMat(c,r);
testCase.assertEqual(full(TM + pi), T + pi);
testCase.assertEqual(full(TM - pi), T - pi);
testCase.assertEqual(full(pi + TM), pi + T);
testCase.assertEqual(full(pi - TM), pi - T);

s = randn(1,1);
TM = s + TM;
TM = TM - s;
testCase.assertEqual(full(TM), T, 'AbsTol', 10*eps, 'RelTol', 10*eps);

end

function test_uminus(testCase)
T = ToepMat([],[]);
testCase.assertTrue(isempty(full(-T)));

T = ToepMat([0,0,0,0], [0,0,0,0]);
testCase.assertEqual(full(-T), zeros(4));

T = ToepMat(4,4);
testCase.assertEqual(full(-T), -4);

[c,r,T] = random_toeplitz(7,7);
TM = ToepMat(c,r);
testCase.assertEqual(full(-TM), -T);
end

function test_uplus(testCase)
% Unary + does not alter the data in any way
[c,r, T] = random_toeplitz(5, 5);
TM = ToepMat(c,r);
testCase.assertEqual(full(+TM), T);
end

function test_add_dense(testCase)

% [] is a Toeplitz matrix, so we stay in this class
T = ToepMat([],[]);
B = [] + T;
testCase.assertEqual(class(B), 'ToepMat');
testCase.assertTrue(isempty(full(B)));

T = toepeye(3);
testCase.assertError( @() T + rand(2), 'tlzstein:InconsistentInput');
testCase.assertError( @() rand(2) - T, 'tlzstein:InconsistentInput');

[c,r,T] = random_toeplitz(9,9);
TM = ToepMat(c,r);
A = randn(9);
testCase.assertEqual(class(TM + A), 'double');
testCase.assertEqual(class(TM - A), 'double');
testCase.assertEqual(class(A + TM), 'double');
testCase.assertEqual(class(A - TM), 'double');
testCase.assertEqual(A + TM, A + T);
testCase.assertEqual(A - TM, A - T);
testCase.assertEqual(TM + A, T + A);
testCase.assertEqual(TM - A, T - A);

end

function test_add_dense_toeplitz(testCase)
% Magic op, where the other op is dense toeplitz.  We notice and convert to
% Toeplitz.
T = toepeye(3);
T = ones(3) - T; % Result is Toeplitz
testCase.assertEqual(class(T), 'ToepMat');
testCase.assertEqual(full(T), toeplitz([0,1,1]));

[c,r,T] = random_toeplitz(7,7);
TM = ToepMat(c,r);
Z = zeros(7);
testCase.assertEqual(class(TM + Z), 'ToepMat');
testCase.assertEqual(full(TM + Z), T);
testCase.assertEqual(class(TM - Z), 'ToepMat');
testCase.assertEqual(full(TM - Z), T);
testCase.assertEqual(class(Z + TM), 'ToepMat');
testCase.assertEqual(full(Z + TM), T);
testCase.assertEqual(class(Z - TM), 'ToepMat');
testCase.assertEqual(full(Z - TM), -T);

[c,r,T] = random_toeplitz(7,7);
TM = ToepMat(c,r);
E = -3*eye(7);
testCase.assertEqual(class(TM + E), 'ToepMat');
testCase.assertEqual(full(TM + E), T + E);
testCase.assertEqual(class(TM - E), 'ToepMat');
testCase.assertEqual(full(TM - E), T - E);
testCase.assertEqual(class(E + TM), 'ToepMat');
testCase.assertEqual(full(E + TM), E + T);
testCase.assertEqual(class(E - TM), 'ToepMat');
testCase.assertEqual(full(E - TM), E - T);

[c,r,T] = random_toeplitz(7,7);
TM = ToepMat(c,r);
cc = rand(7,1);
rr = rand(7,1);
cc(1) = rr(1);
A = toeplitz(cc,rr);
testCase.assertEqual(class(TM + A), 'ToepMat');
testCase.assertEqual(full(TM + A), T + A);
testCase.assertEqual(class(TM - A), 'ToepMat');
testCase.assertEqual(full(TM - A), T - A);
testCase.assertEqual(class(A + TM), 'ToepMat');
testCase.assertEqual(full(A + TM), A + T);
testCase.assertEqual(class(A - TM), 'ToepMat');
testCase.assertEqual(full(A - TM), A - T);

end

function test_add_toepmat(testCase)
T = ToepMat([], []) + ToepMat([], []);
testCase.assertEqual(class(T), 'ToepMat');
testCase.assertTrue(isempty(full(T)));

T1 = ToepMat(3,3);
T2 = ToepMat(-3,-3);
testCase.assertEqual(class(T1+T2), 'ToepMat');
testCase.assertEqual(class(T1-T2), 'ToepMat');
testCase.assertEqual(class(T2+T1), 'ToepMat');
testCase.assertEqual(class(T2-T1), 'ToepMat');
testCase.assertEqual(full(T1 + T2), 0);
testCase.assertEqual(full(T1 - T2), 6);
testCase.assertEqual(full(T2 + T1), 0);
testCase.assertEqual(full(T2 - T1), -6);

[c1,r1,T1] = random_toeplitz(6,6);
[c2,r2,T2] = random_toeplitz(6,6);
TM1 = ToepMat(c1,r1);
TM2 = ToepMat(c2,r2);
testCase.assertEqual(class(TM1+TM2), 'ToepMat');
testCase.assertEqual(class(TM1-TM2), 'ToepMat');
testCase.assertEqual(class(TM2+TM1), 'ToepMat');
testCase.assertEqual(class(TM2-TM1), 'ToepMat');

testCase.assertEqual(full(TM1 + TM2), T1+T2);
testCase.assertEqual(full(TM1 - TM2), T1-T2);
testCase.assertEqual(full(TM2 + TM1), T2+T1);
testCase.assertEqual(full(TM2 - TM1), T2-T1);

TM1 = ToepMat([1,2,3], [1,5,6]);
TM2 = ToepMat([6,5], [6,0]);
testCase.assertError( @() TM1 + TM2, 'tlzstein:InconsistentInput');


end

function test_add_tlmat(testCase)
% Result is always TLMat.

n = 11;
[c,r,T] = random_toeplitz(n,n);
TM = ToepMat(c,r);
TL = TLMat(rand(n,3), rand(n,3));

testCase.assertEqual(class(TM + TL), 'TLMat');
testCase.assertEqual(class(TM - TL), 'TLMat');
testCase.assertEqual(class(TL + TM), 'TLMat');
testCase.assertEqual(class(TL - TM), 'TLMat');

% May fail with prob 0
testCase.assertEqual(drank(TM + TL), 5);
testCase.assertEqual(drank(TM - TL), 5);
testCase.assertEqual(drank(TL + TM), 5);
testCase.assertEqual(drank(TL - TM), 5);

testCase.assertEqual(full(TM + TL), T + full(TL), ...
    'RelTol', 100*eps, 'AbsTol', 100*eps);
testCase.assertEqual(full(TM - TL), T - full(TL), ...
    'RelTol', 100*eps, 'AbsTol', 100*eps);
testCase.assertEqual(full(TL + TM), full(TL) + T, ...
    'RelTol', 100*eps, 'AbsTol', 100*eps);
testCase.assertEqual(full(TL - TM), full(TL) - T, ...
    'RelTol', 100*eps, 'AbsTol', 100*eps);

end


function test_mtimes_scalar(testCase)

TM = ToepMat([], []);
s = 0;
testCase.assertEqual(class(TM * s), 'ToepMat');
testCase.assertEqual(class(s * TM), 'ToepMat');
testCase.assertTrue(isempty(full(TM * s)));
testCase.assertTrue(isempty(full(s * TM)));

s = 4 - 2i;
testCase.assertEqual(class(TM * s), 'ToepMat');
testCase.assertEqual(class(s * TM), 'ToepMat');
testCase.assertTrue(isempty(full(TM * s)));
testCase.assertTrue(isempty(full(s * TM)));

t = -1 + sqrt(2)*1i;
TM = ToepMat(t,t);
s = 0;
testCase.assertEqual(class(TM * s), 'ToepMat');
testCase.assertEqual(class(s * TM), 'ToepMat');
testCase.assertEqual(full(TM * s), 0);
testCase.assertEqual(full(s * TM), 0);

s = 4 - 2i;
testCase.assertEqual(class(TM * s), 'ToepMat');
testCase.assertEqual(class(s * TM), 'ToepMat');
testCase.assertEqual(full(TM * s), t*s);
testCase.assertEqual(full(s * TM), s*t);

TM = toepeye(8);
testCase.assertEqual(class(TM * s), 'ToepMat');
testCase.assertEqual(class(s * TM), 'ToepMat');
testCase.assertEqual(full(TM * s), s*eye(8));
testCase.assertEqual(full(s * TM), s*eye(8));

[c,r,T] = random_toeplitz(12,12);
TM = ToepMat(c,r);
s = exp(12/5 * pi * 1i);
testCase.assertEqual(class(TM * s), 'ToepMat');
testCase.assertEqual(class(s * TM), 'ToepMat');
testCase.assertEqual(full(TM * s), T*s);
testCase.assertEqual(full(s * TM), s*T);

end

function test_mtimes_tritoep(testCase)
% Magic case where both operands are triangular, product stays in Toeplitz
% space.
end

function test_mtimes_toep(testCase)
% Both operands are ToepMat
end

function test_mtimes_tl(testCase)
% One op is a TLMat
end