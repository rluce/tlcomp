function tests = ToepMatTest

tests = functiontests(localfunctions);

end


function test_constructor_empty(testCase)
T = ToepMat([], []);
testCase.assertTrue(isempty(T.c));
testCase.assertTrue(isempty(T.r));
end

function test_construct_scalar(testCase)
T = ToepMat(-1, -1);
testCase.assertEqual(T.c, -1);
testCase.assertEqual(T.r, -1);

T = ToepMat(3);
testCase.assertEqual(T.c, 3);
testCase.assertEqual(T.r, 3);

% Conflicting information for diagonal element.  We align with the Matlab
% built-in 'toeplitz' and let the column information win.
c = 6;
r = 1;
id = 'tlcomp:InconsistentInput';
testCase.assertWarning( @() ToepMat(c,r), id);
warning('off', id);
T = ToepMat(c,r);
warning('on', id);
testCase.assertEqual(T.c, 6);
testCase.assertEqual(T.r, 6);

end


function test_constructor_two_vec(testCase)
[c,r] = random_toeplitz(7,7);

T = ToepMat(c,r);
testCase.assertEqual(T.c, c(:));
testCase.assertEqual(T.r, r(:));

% Conflicting information for diagonal element.  We align with the Matlab
% built-in 'toeplitz' and let the column information win.
c = rand(9,1);
r = rand(9,1);
id = 'tlcomp:InconsistentInput';
testCase.assertWarning( @() ToepMat(c,r), id);
warning('off', id);
T = ToepMat(c,r);
warning('on', id);
testCase.assertEqual(T.c, c(:));
testCase.assertEqual(T.r, [c(1); r(2:end)]);


% For the moment, we only support square matrices
c = ones(4,1);
r = ones(3,1);
testCase.assertError( @() ToepMat(c,r), 'tlcomp:InconsistentInput');

% Only vectors allowed
c = ones(4,2);
r = ones(4,2);
testCase.assertError( @() ToepMat(c,r), 'tlcomp:InconsistentInput');

end

function test_construct_one_vec(testCase)
% Construct a Hermitian Toeplitz matrix
n = 7;
r = randn(n,1) + 1i * randn(n,1);
r(1) = abs(r(1));
T = ToepMat(r);
testCase.assertEqual(T.r, r(:));
testCase.assertEqual(T.c, conj(r(:)));
testCase.assertEqual(full(T), toeplitz(r)); % Align with Matlab built-in

% Diagonal not real -> Hermitian only on the off-diagonals
r = randn(n,1) + 1i * randn(n,1);
T = ToepMat(r);
testCase.assertEqual(T.r, r(:));
testCase.assertEqual(T.c, [r(1); conj(r(2:end))]);
testCase.assertEqual(full(T), toeplitz(r)); % Align with Matlab built-in


% Only vectors allowed
r = ones(4,2);
testCase.assertError( @() ToepMat(r), 'tlcomp:InconsistentInput');

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
testCase.assertEqual(size(T, 1), 3);
testCase.assertEqual(size(T, 2), 3);


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
testCase.assertError( @() T + rand(2), 'tlcomp:InconsistentInput');
testCase.assertError( @() rand(2) - T, 'tlcomp:InconsistentInput');

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
testCase.assertError( @() TM1 + TM2, 'tlcomp:InconsistentInput');


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

function test_mtimes_toep(testCase)
% Both operands are ToepMat, the result is TLMat
% TODO This could be done a bit better:  if both a triangular, the result
% is Toeplitz again, for example.  For now, the result will always be TL.

TM1 = ToepMat([1,1], [1,2]);
TM2 = ToepMat([1,1,1], [1,2,2]);
testCase.assertError( @() TM1 * TM2, 'tlcomp:InconsistentInput');


TM1 = ToepMat([], []);
TM2 = TM1;
testCase.assertEqual(class(TM1 * TM2), 'TLMat');
testCase.assertTrue(isempty(full(TM1 *TM2)));

TM1 = ToepMat(4,4);
TM2 = ToepMat(.5, .5);
testCase.assertEqual(class(TM1 * TM2), 'TLMat');
testCase.assertEqual(full(TM1 * TM2), 2);
testCase.assertEqual(class(TM2 * TM1), 'TLMat');
testCase.assertEqual(full(TM2 * TM1), 2);

[c, r, T] = random_toeplitz(8,8);
TM1 = ToepMat(c,r);
TM2 = toepeye(8);
nfact = 8*norm(full(TM1)) * norm(full(TM2));
testCase.assertEqual(class(TM1 * TM2), 'TLMat');
testCase.assertEqual(full(TM1 * TM2), T, ...
    'RelTol', nfact*eps, 'AbsTol', nfact*eps);
testCase.assertEqual(class(TM2 * TM1), 'TLMat');
testCase.assertEqual(full(TM2 * TM1), T, ...
    'RelTol', nfact*eps, 'AbsTol', nfact*eps);
testCase.assertEqual(drank(TM1*TM2), 2);
testCase.assertEqual(drank(TM2*TM1), 2);

[c, r, T] = random_toeplitz(8,8);
TM1 = ToepMat(c,r);
TM2 = ToepMat(ones(8,1), ones(8,1));
testCase.assertEqual(class(TM1 * TM2), 'TLMat');
testCase.assertEqual(full(TM1 * TM2), T * ones(8), ...
    'RelTol', 500*eps, 'AbsTol', 500*eps);
testCase.assertEqual(class(TM2 * TM1), 'TLMat');
testCase.assertEqual(full(TM2 * TM1), ones(8) * T, ...
    'RelTol', 500*eps, 'AbsTol', 500*eps);
testCase.assertEqual(drank(TM1*TM2), 2);
testCase.assertEqual(drank(TM2*TM1), 1);

[c1, r1, T1] = random_toeplitz(9,9);
TM1 = ToepMat(c1,r1);
[c2, r2, T2] = random_toeplitz(9,9);
TM2 = ToepMat(c2,r2);
testCase.assertEqual(class(TM1 * TM2), 'TLMat');
testCase.assertEqual(full(TM1 * TM2), T1 * T2, ...
    'RelTol', 500*eps, 'AbsTol', 500*eps);
testCase.assertEqual(class(TM2 * TM1), 'TLMat');
testCase.assertEqual(full(TM2 * TM1), T2 * T1, ...
    'RelTol', 500*eps, 'AbsTol', 500*eps);
testCase.assertEqual(drank(TM1*TM2), 4);
testCase.assertEqual(drank(TM2*TM1), 4);



end

function test_mtimes_tl(testCase)
% One op is a TLMat

TM = ToepMat([1,1], [1,2]);
TL = TLMat([1,1,1], [1,2,2]);
testCase.assertError( @() TM * TL, 'tlcomp:InconsistentInput');

TM = ToepMat([], []);
TL = TLMat([]);
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertTrue(isempty(full(TM *TL)));

TM = ToepMat(4,4);
TL = TLMat(.5, .5);
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertEqual(full(TM * TL), 2);
testCase.assertEqual(class(TL * TM), 'TLMat');
testCase.assertEqual(full(TL * TM), 2);

[c, r, T] = random_toeplitz(8,8);
TM = ToepMat(c,r);
TL = tleye(8);
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertEqual(full(TM * TL), T, 'AbsTol', 100*eps, 'RelTol', 100*eps);
testCase.assertEqual(class(TL * TM), 'TLMat');
testCase.assertEqual(full(TL * TM), T, 'AbsTol', 50*eps, 'RelTol', 50*eps);
testCase.assertEqual(drank(TM*TL), 2);
testCase.assertEqual(drank(TL*TM), 2);

[c, r, T] = random_toeplitz(8,8);
TM = ToepMat(c,r);
TL = TLMat(ones(8,1), ones(8,1));
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertEqual(full(TM * TL), T * ones(8), ...
    'AbsTol', 500*eps, 'RelTol', 500*eps);
testCase.assertEqual(class(TL * TM), 'TLMat');
testCase.assertEqual(full(TL * TM), ones(8) * T, ...
    'AbsTol', 500*eps, 'RelTol', 500*eps);
testCase.assertEqual(drank(TM*TL), 2);
testCase.assertEqual(drank(TL*TM), 1);

[c1, r1, T1] = random_toeplitz(9,9);
TM = ToepMat(c1,r1);
[c2, r2, T2] = random_toeplitz(9,9);
TL = TLMat(c2,r2);
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertEqual(full(TM * TL), T1 * T2, 'AbsTol', 500*eps, 'RelTol', 500*eps);
testCase.assertEqual(class(TL * TM), 'TLMat');
testCase.assertEqual(full(TL * TM), T2 * T1, 'AbsTol', 500*eps, 'RelTol', 500*eps);
testCase.assertEqual(drank(TM*TL), 4);
testCase.assertEqual(drank(TL*TM), 4);

[c, r, T] = random_toeplitz(12,12);
TM = ToepMat(c,r);
TL = TLMat(rand(12,4),rand(12,4));

nfact = 8 * norm(full(TM)) * norm(full(TL));
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertEqual(full(TM * TL), T * full(TL), ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(class(TL * TM), 'TLMat');
testCase.assertEqual(full(TL * TM), full(TL) * T, ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(drank(TM*TL), 6);
testCase.assertEqual(drank(TL*TM), 6);


end

function test_mtimes_double_matrix(testCase)
% One operand is a dense matrix, result is dense
% TODO one could do better and detect Toeplitz structure.  For now just
% return dense.

TM = ToepMat([1,1], [1,2]);
A = randn(3);
testCase.assertError( @() TM * A, 'tlcomp:InconsistentInput');
testCase.assertError( @() A * TM, 'tlcomp:InconsistentInput');

TM = ToepMat([],[]);
A = [];
testCase.assertEqual(class(TM * A), 'double');
testCase.assertEqual(TM * A, []);
testCase.assertEqual(class(A * TM), 'double');
testCase.assertEqual(A * TM, []);

% This should resolve through scalar multiplication
TM = ToepMat(1,1);
A = 1i;
testCase.assertEqual(class(TM * A), 'ToepMat');
testCase.assertEqual(full(TM * A), 1i);
testCase.assertEqual(class(A * TM), 'ToepMat');
testCase.assertEqual(full(A * TM), 1i);

[c,r,T] = random_toeplitz(7,7);
TM = ToepMat(c,r);
A = eye(7);
testCase.assertEqual(class(TM * A), 'double');
testCase.assertEqual(TM * A, T, 'AbsTol', 10*eps, 'RelTol', 10*eps);
testCase.assertEqual(class(A * TM), 'double', 'AbsTol', 10*eps, 'RelTol', 10*eps);
testCase.assertEqual(A * TM, T, 'AbsTol', 10*eps, 'RelTol', 10*eps);

[c,r,T] = random_toeplitz(8,8);
TM = ToepMat(c,r);
A = ones(8);
testCase.assertEqual(class(TM * A), 'double');
testCase.assertEqual(TM * A, sum(T,2) * ones(1,8), ...
    'AbsTol', 10*eps, 'RelTol', 10*eps);
testCase.assertEqual(class(A * TM), 'double');
testCase.assertEqual(A * TM, ones(8,1) * sum(T,1), ...
    'AbsTol', 10*eps, 'RelTol', 10*eps);

[c,r,T] = random_toeplitz(8,8);
TM = ToepMat(c,r);
A = rand(8);

nfact = norm(A) * norm(full(TM));
testCase.assertEqual(class(TM * A), 'double');
testCase.assertEqual(TM * A, T * A, 'RelTol', 4*nfact*eps, 'RelTol', 4*nfact*eps);
testCase.assertEqual(class(A * TM), 'double');
testCase.assertEqual(A * TM, A * T, 'AbsTol', 4*nfact*eps, 'RelTol', 4*nfact*eps);
end

function test_mtimes_double_vector(testCase)
% Important special case: Toeplitz times vector

[c,r,T] = random_toeplitz(8,8);
TM = ToepMat(c,r);
x = zeros(8,1);
x(3) = 1;
testCase.assertEqual(class(TM * x), 'double');
testCase.assertEqual(TM * x, T(:,3), 'AbsTol', 10*eps, 'RelTol', 10*eps);
testCase.assertEqual(class(x' * TM), 'double');
testCase.assertEqual(x' * TM, T(3,:), 'AbsTol', 10*eps, 'RelTol', 10*eps);

[c,r,T] = random_toeplitz(8,8);
TM = ToepMat(c,r);
x = ones(8,1);
testCase.assertEqual(class(TM * x), 'double');
testCase.assertEqual(TM * x, sum(T,2), 'RelTol', 10*eps, 'AbsTol', 10*eps);
testCase.assertEqual(class(x' * TM), 'double');
testCase.assertEqual(x' * TM, sum(T,1), 'AbsTol', 10*eps, 'RelTol', 10*eps);

[c,r,T] = random_toeplitz(12,12);
TM = ToepMat(c,r);
x = rand(12,1);
nfact = norm(full(TM)) * norm(x);
testCase.assertEqual(class(TM * x), 'double');
testCase.assertEqual(TM * x, T*x, 'AbsTol', 10*eps, 'RelTol', 10*eps);
testCase.assertEqual(class(x' * TM), 'double');
testCase.assertEqual(x' * TM, x'*T, 'AbsTol', 4*nfact*eps, 'RelTol', 4*nfact*eps);
end

function test_transpose(testCase)

TM = ToepMat([], []);
testCase.assertEqual(class(TM.'), 'ToepMat');
testCase.assertEqual(full(TM.'), []);

TM = ToepMat(2,2);
testCase.assertEqual(class(TM.'), 'ToepMat');
testCase.assertEqual(full(TM.'), 2);

TM = ToepMat(2i,2i);
testCase.assertEqual(class(TM.'), 'ToepMat');
testCase.assertEqual(full(TM.'), 2i);

[c,r,T] = random_toeplitz(9,9);
TM = ToepMat(c,r);
testCase.assertEqual(class(TM.'), 'ToepMat');
testCase.assertEqual(full(TM.'), T.');
testCase.assertEqual(full((TM.').'), T);
end

function test_ctranspose(testCase)

TM = ToepMat([], []);
testCase.assertEqual(class(TM'), 'ToepMat');
testCase.assertEqual(full(TM'), []);

TM = ToepMat(2,2);
testCase.assertEqual(class(TM'), 'ToepMat');
testCase.assertEqual(full(TM'), 2);

TM = ToepMat(2i,2i);
testCase.assertEqual(class(TM'), 'ToepMat');
testCase.assertEqual(full(TM'), -2i);

[c,r,T] = random_toeplitz(9,9);
TM = ToepMat(c,r);
testCase.assertEqual(class(TM'), 'ToepMat');
testCase.assertEqual(full(TM'), T');
testCase.assertEqual(full((TM')'), T);
end

function test_polyvalm(testCase)

p = [1, 1, 1];
TM = 1i * toepeye(9);
T = 1i * eye(9);

P = polyvalm(p, TM);
P_true = polyvalm(p, T);
testCase.assertEqual(class(P), 'TLMat');
testCase.assertEqual(full(P), P_true, 'AbsTol', 2*eps, 'RelTol', 2*eps);

p = randn(1, 2) + 1i * randn(1,2);
[c,r,T] = random_toeplitz(12,12);
TM = ToepMat(c,r);
P = polyvalm(p, TM);
P_true = polyvalm(p, T);
testCase.assertEqual(class(P), 'TLMat');
testCase.assertEqual(full(P), P_true, 'AbsTol', 100*eps, 'RelTol', 100*eps);

end

function test_mldivide_double(testCase)
% Toeplitz \ (double matrix/vector)

TM = 1i * toepeye(1);
testCase.assertEqual(TM\0, 0);
testCase.assertEqual(TM\1, -1i);
testCase.assertEqual(TM\1i, 1);

TM = toepeye(7);
B = randn(7,5) + 1i * randn(7,5);
X = TM \ B;
testCase.assertEqual(X, B, 'AbsTol', 10*eps, 'RelTol', 10*eps);

[c,r,T] = random_toeplitz(12,12);
x = rand(12,1);
b = T*x;
TM = ToepMat(c,r);
testCase.assertEqual(TM\b, complex(x), 'AbsTol', 50*eps, 'RelTol', 50*eps);

[c,r,T] = random_toeplitz(12,12);
x = rand(12,5) + 1i * rand(12,5);
b = T*x;
TM = ToepMat(c,r);
testCase.assertEqual(TM\b, x, 'AbsTol', 100*eps, 'RelTol', 100*eps);

[c,r,T] = random_toeplitz(12,12);
c = real(c); r = real(r); T = real(T);
x = rand(12,5);
b = T*x;
TM = ToepMat(c,r);
testCase.assertTrue(isreal(TM\b));
nfact = 8 * cond(full(TM)) * norm(b);
testCase.assertEqual(TM\b, x, 'AbsTol', nfact*eps, 'RelTol', nfact*eps);

TM = toepeye(8);
testCase.assertError( @() TM \ 1.0, 'tlcomp:InconsistentInput');
end

function test_mldivide_toepmat(testCase)

TM1 = toepeye(9);
TM2 = toepeye(9);

TM = TM1 \ TM2;
% We could do better and detect truely Toeplitz matrices, but for now we
% always promote to TLMat
testCase.assertTrue(isa(TM, 'TLMat'));
testCase.assertEqual(full(TM), eye(9), 'AbsTol', eps, 'RelTol', eps);

% Cayley transform
e2 = zeros(16,1);
e2(2) = 1;
TM = ToepMat(-e2, e2); % Skew-symm
T = full(TM);
R = (TM - toepeye(16)) \ (TM + toepeye(16));
testCase.assertEqual(full(R), (T - eye(16)) \ (T + eye(16)), ...
    'AbsTol', 50*eps, 'RelTol', 50*eps);

[c1, r1, T1] = random_toeplitz(7,7);
[c2, r2, T2] = random_toeplitz(7,7);
TM1 = ToepMat(c1, r1);
TM2 = ToepMat(c2, r2);
TM = TM1 \ TM2;

nfact = 4*cond(full(TM1)) * norm(full(TM2));
testCase.assertTrue(isa(TM, 'TLMat'));
testCase.assertEqual(full(TM), T1 \T2, 'AbsTol', nfact*eps, 'RelTol', nfact*eps);

end

function test_ratvalm(testCase)

% A fancy way to write the identity mapping
p = [1, 0];
q = 1;
[c, r, T] = random_toeplitz(13,13);
TM = ToepMat(c,r);
rT = TM.ratevalm(p,q);
testCase.assertTrue(isa(rT, 'TLMat'));
testCase.assertEqual(full(rT), T, 'AbsTol', 50*eps, 'RelTol', 50*eps);

% Cayley transform
e2 = zeros(16,1);
e2(2) = 1;
TM = ToepMat(-e2, e2); % Skew-symm
T = full(TM);
p = [1, 1];
q = [1, -1];
rT = TM.ratevalm(p, q);
testCase.assertEqual(full(rT), (T - eye(16)) \ (T + eye(16)), ...
    'AbsTol', 50*eps, 'RelTol', 50*eps);

% Random rational function
p = randn(1,4) + 1i * randn(1,4);
q = randn(1,4) + 1i * randn(1,4);
[c,r,T] = random_toeplitz(9,9);
TM = ToepMat(c,r);
rT = TM.ratevalm(p,q);
nfact = cond(polyvalm(q, T)) * norm(polyvalm(p, T));
testCase.assertEqual(full(rT), polyvalm(q, T) \ polyvalm(p, T), ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
end


function test_inv(testCase)

TM = toepeye(3);
TMinv = inv(TM);
testCase.assertEqual(class(TMinv), 'TLMat');
testCase.assertEqual(full(TMinv), eye(3), 'AbsTol', 10*eps, 'RelTol', 10*eps);

[c,r,T] = random_toeplitz(9,9);
TM = ToepMat(c,r);
testCase.assertEqual(full(inv(TM)), inv(T), 'AbsTol', ...
    1e4*eps, 'RelTol', 1e4*eps);

end

function test_norm(testCase)

[c, r, T] = random_toeplitz(7,7);
TM = ToepMat(c,r);
testCase.assertEqual(norm(TM, 1), norm(T, 1), 'RelTol', 4*eps, 'AbsTol', 4*eps);
testCase.assertEqual(norm(TM, 'inf'), norm(T, 'inf'), 'RelTol', 4*eps, 'AbsTol', 4*eps);
testCase.assertEqual(norm(TM, inf), norm(T, inf), 'RelTol', 4*eps, 'AbsTol', 4*eps);
testCase.assertEqual(norm(TM, 'fro'), norm(T, 'fro'), 'RelTol', 4*eps, 'AbsTol', 4*eps);

testCase.assertWarning(@() norm(TM), 'tlcomp:Unsupported');
testCase.assertWarning(@() norm(TM, 2), 'tlcomp:Unsupported');
testCase.assertError(@() norm(TM, 3), 'tlcomp:InconsistentInput');

warning('OFF', 'tlcomp:Unsupported');
testCase.assertEqual(norm(TM), norm(T));
testCase.assertEqual(norm(TM, 2), norm(T, 2));
warning('ON', 'tlcomp:Unsupported');

end

function test_mrdivide(testCase)

TM = toepeye(1);
A = TM / 1;
testCase.assertTrue(isa(A, 'ToepMat'));
testCase.assertEqual(full(A), 1);

TM = toepeye(6);
A = TM / 4;
testCase.assertTrue(isa(A, 'ToepMat'));
testCase.assertEqual(full(A), .25 * eye(6));

[c, r, T] = random_toeplitz(9,9);
TM = ToepMat(c,r);
testCase.assertEqual(full(TM/1), T/1);
testCase.assertEqual(full(TM/2), T/2);
testCase.assertEqual(full(TM/5), T/5, 'AbsTol', 10*eps, 'RelTol', 10*eps);
testCase.assertEqual(full(TM/.2), T/.2, 'AbsTol', 10*eps, 'RelTol', 10*eps);
end

function test_mpower_monomial(testCase)
% Note: We would need to work in ToepMat to detect whether a T^s is in fact
% a Toeplitz matrix.  We take an easy design choice and guaratee ToepMat
% for s=0,1 and resort to TLMat for all s > 1.
TM = ToepMat(0);
testCase.assertEqual(full(TM^0), 1);
testCase.assertTrue(isa(TM^0, 'ToepMat'));
testCase.assertEqual(full(TM^1), 0);
testCase.assertTrue(isa(TM^1, 'ToepMat'));
testCase.assertEqual(full(TM^2), 0);
testCase.assertTrue(isa(TM^2, 'TLMat'));

TM = toepeye(1);
for s = 0:4
    testCase.assertEqual(full(TM^s), 1.0);
end

TM = toepeye(3);
for s = 0:4
    testCase.assertEqual(full(TM^s), eye(3))
end
testCase.assertTrue(isa(TM^0, 'ToepMat'));
testCase.assertTrue(isa(TM^1, 'ToepMat'));
testCase.assertTrue(isa(TM^2, 'TLMat'));

[c,r,T] = random_toeplitz(8,8);
TM = ToepMat(c,r);
for s = 0:4
    Tpow_true = T^s;
    Tpow = TM^s;
    nT = norm(Tpow_true);
    testCase.assertEqual(full(Tpow), Tpow_true, ...
        'AbsTol', 4*nT*eps, 'RelTol', 4*nT*eps);
end
testCase.assertTrue(isa(TM^0, 'ToepMat'));
testCase.assertTrue(isa(TM^1, 'ToepMat'));
testCase.assertTrue(isa(TM^2, 'TLMat'));


end


function test_det(testCase)

T = ToepMat(zeros(5,1));
testCase.assertEqual(det(T), 0.0);

T = pi * toepeye(3);
testCase.assertEqual(det(T), pi^3);

c = randn(5,1);
r = randn(5,1);
r(1) = c(1);
T = ToepMat(c, r);
testCase.assertEqual(det(T), det(full(T)), 'RelTol', 256 * eps);
end