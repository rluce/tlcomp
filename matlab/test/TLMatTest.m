function tests = TLMatTest

tests = functiontests(localfunctions);

end


function test_construct_empty(testCase)
TL = TLMat([]);
testCase.assertTrue(isempty(TL.G));
testCase.assertTrue(isempty(TL.B));

TL = TLMat([],[]);
testCase.assertTrue(isempty(TL.G));
testCase.assertTrue(isempty(TL.B));

TL = TLMat([],[], 'GB');
testCase.assertTrue(isempty(TL.G));
testCase.assertTrue(isempty(TL.B));

end


function test_construct_scalar(testCase)
% Recall reconstruction formula for 1-by-1 case: T(G,B) = 0.5 * G * B'

A = 6;
TL = TLMat(A);
testCase.assertEqual(0.5 * TL.G * TL.B', A);

% Generators for A = 6
G = 12;
B = 1;
TL = TLMat(G,B, 'GB');
testCase.assertEqual(0.5 * TL.G * TL.B', A);

% Non-minimal 1-by-1 generators for A = 6
G = [2, 2];
B = [3, 3];
TL = TLMat(G,B, 'GB');
testCase.assertEqual(0.5 * TL.G * TL.B', A);

c = 6;
r = 6;
TL = TLMat(c,r);
testCase.assertEqual(0.5 * TL.G * TL.B', A);

c = 6;
r = 1;
testCase.assertError( @() TLMat(c,r), 'tlcomp:InconsistentInput');


end


function test_construct_eye(testCase)
n = 8;

E = eye(8);
e1 = E(:,1);
en = E(:,end);
z = zeros(n,1);

TL = TLMat(E);
testCase.assertEqual(TL.G * TL.B', 2*e1*en');

TL = TLMat(e1,e1);
testCase.assertEqual(TL.G * TL.B', 2*e1*en');

TL = TLMat(2*e1, en, 'GB');
testCase.assertEqual(TL.G * TL.B', 2*e1*en');

TL = TLMat([2*e1,z], [en,z], 'GB');
testCase.assertEqual(TL.G * TL.B', 2*e1*en');

TL = TLMat([2*e1,z], [en,z]);
testCase.assertEqual(TL.G * TL.B', 2*e1*en');

end

function test_construct_random(testCase)
n = 9;
[c,r,T] = random_toeplitz(n,n);

TL = TLMat(T);
T2 = toeplkreconstruct(TL.G, TL.B);
testCase.assertEqual(T2, T, 'AbsTol', 100*eps, 'RelTol', 100*eps);

TL = TLMat(c,r);
T2 = toeplkreconstruct(TL.G, TL.B);
testCase.assertEqual(T2, T, 'AbsTol', 100*eps, 'RelTol', 100*eps);

[G, B] = toepgen(c,r);
TL = TLMat(G,B);
T2 = toeplkreconstruct(TL.G, TL.B);
testCase.assertEqual(T2, T, 'AbsTol', 100*eps, 'RelTol', 100*eps);

GG = [2 * G, -G];
BB = [B, B];
TL = TLMat(GG, BB);
T2 = toeplkreconstruct(TL.G, TL.B);
testCase.assertEqual(T2, T, 'AbsTol', 100*eps,  'RelTol', 100*eps);


end

function test_construct_badinput(testCase)

testCase.assertError( @() TLMat(rand(4,2), rand(4,3)), ...
    'tlcomp:InconsistentInput');

end


function test_size(testCase)

TL = TLMat([]);
testCase.assertEqual(size(TL), [0,0]);

TL = TLMat(1);
testCase.assertEqual(size(TL), [1,1]);

TL = TLMat([1,2,5], [1,2,5]);
testCase.assertEqual(size(TL), [3,3]);

testCase.assertEqual(size(TL, 1), 3);
testCase.assertEqual(size(TL, 2), 3);

end


function test_drank(testCase)

TL = TLMat([]);
testCase.assertEqual(TL.drank, 0);

TL = TLMat(eye(11));
testCase.assertEqual(TL.drank, 1);

[c,r] = random_toeplitz(9,9);
TL = TLMat(c,r);
testCase.assertEqual(TL.drank, 2);

TL = TLMat(rand(9));
testCase.assertEqual(TL.drank, 9);

end

function test_full(testCase)

E = eye(13);
TL = TLMat(E(:,1), E(:,1));
testCase.assertEqual(full(TL), E);

[c, r, T] = random_toeplitz(5,5);
TL = TLMat(c,r);
testCase.assertEqual(full(TL), T, 'RelTol', 100 * eps);

end


function test_plus_scalar(testCase)

TL = TLMat([0,0], [0,0]);
TLp1 = TL + 1;
testCase.assertEqual(full(TLp1), ones(2), 'AbsTol', 5 * eps);

[c, r] = random_toeplitz(7,7);
TL = TLMat(c,r);

TLp0 = 0 + TL + 0;
testCase.assertEqual(full(TLp0), full(TL));

TLp1 = TL + 1;
testCase.assertEqual(full(TLp1), full(TL) + 1, ...
    'RelTol', 100*eps, 'AbsTol', 100*eps);

TLmpi = TL - pi;
testCase.assertEqual(full(TLmpi), full(TL) - pi, 'RelTol', 100*eps);

n = 9;
e = ones(n,1);
TL = TLMat(e,e);
TL = TL + 1;
testCase.assertEqual(full(TL), 2*ones(n), 'RelTol', 10*eps);
TL = -2 + TL;
testCase.assertEqual(full(TL), zeros(n), 'AbsTol', 10*eps);

end


function test_plus_dense_matrix(testCase)

TL = TLMat([]);
A = [];

% Adding empty matrix to empty TL matrix stays empty TL matrix
B = A + TL;
testCase.assertEqual(class(B), 'TLMat');
testCase.assertTrue(isempty(full(B)));

% This is in fact scalar addition
TL = tleye(1);
A = eye(1);
B = TL + A;
testCase.assertEqual(class(B), 'TLMat');
testCase.assertEqual(full(B), 2);

TL = tleye(2);
E2 = eye(2);
S = E2 - TL;
testCase.assertEqual(full(S), zeros(2));

n = 7;
[c,r] = random_toeplitz(n,n);
TL = TLMat(c,r);
A = randn(n,n);
B = TL + A;
testCase.assertEqual(B, full(TL) + A);
testCase.assertEqual(class(B), 'double');

end

function test_plus_dense_toeplitz(testCase)
% Magic function:  Detect wether other operand is a dense toeplitz matrix

n = 13;
TL = TLMat(1i * rand(n,3) + rand(n,3), rand(n,3));
nTL = norm(full(TL));
B = TL + zeros(n);
testCase.assertEqual(class(B), 'TLMat');
testCase.assertEqual(drank(B), 3);
testCase.assertEqual(full(B), full(TL), 'AbsTol', nTL*eps, 'RelTol', nTL*eps);

sigma = rand(1,1) + 1i * rand(1,1);
B = TL + sigma * eye(n);
testCase.assertEqual(class(B), 'TLMat');
testCase.assertEqual(drank(B), 4); % fails with prob 0
testCase.assertEqual(full(B), full(TL) + sigma * eye(n), ...
    'AbsTol', 100*eps, 'RelTol', 100*eps);

[~, ~, T] = random_toeplitz(n,n);
testCase.assertEqual(class(TL + T), 'TLMat');
testCase.assertEqual(class(TL - T), 'TLMat');
testCase.assertEqual(class(T + TL), 'TLMat');
testCase.assertEqual(class(T - TL), 'TLMat');

testCase.assertEqual(drank(TL + T), 5); % fails with prob 0
testCase.assertEqual(drank(TL - T), 5); % fails with prob 0
testCase.assertEqual(drank(T + TL), 5); % fails with prob 0
testCase.assertEqual(drank(T - TL), 5); % fails with prob 0

nfact = 2 * (norm(full(T)) + norm(full(TL)));
testCase.assertEqual(full(TL + T), full(TL) + T, ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(full(TL - T), full(TL) - T, ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(full(T + TL), T + full(TL), ...
    'Abstol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(full(T - TL), T - full(TL), ...
    'Abstol', nfact*eps, 'RelTol', nfact*eps);

end


function test_plus_tlmatrix(testCase)
TL1 = TLMat([]);
TL2 = TLMat([]);
testCase.assertTrue(isempty(full(TL1 + TL2)));

TL1 = TLMat(1,1);
TL2 = TLMat(-2,-2);
TL = TL1 + TL2;
testCase.assertEqual(drank(TL), 1);
testCase.assertEqual(size(TL), [1,1]);
testCase.assertEqual(full(TL), -1);


[c1, r1] = random_toeplitz(11,11);
[c2, r2] = random_toeplitz(11,11);
TL1 = TLMat(c1, r1);
TL2 = TLMat(c2, r2);
TL = TL1 + TL2;
testCase.assertEqual(drank(TL), 2);
testCase.assertEqual(full(TL), toeplitz(c1,r1) + toeplitz(c2,r2), ...
    'RelTol', 100*eps, 'AbsTol', 100*eps);

end

function test_plus_toepmat(testCase)
% Other operand is a ToepMat
TL = TLMat([], []);
TM = ToepMat([], []);
testCase.assertEqual(class(TL + TM), 'TLMat');
testCase.assertEqual(class(TL - TM), 'TLMat');
testCase.assertEqual(class(TM + TL), 'TLMat');
testCase.assertEqual(class(TM - TL), 'TLMat');


% Important use case: shifts
n = 9;
r = 3;
TL = TLMat(rand(n,r), rand(n,r) + 1i * rand(n,r));
TL_true = full(TL) + diag(ones(n,1));
TL = TL + toepeye(n);

nfact = norm(TL_true);
testCase.assertEqual(class(TL), 'TLMat');
testCase.assertEqual(full(TL), TL_true, 'RelTol', nfact*eps, 'AbsTol', nfact*eps);
testCase.assertEqual(drank(TL), r + 1); % Fails with prob. 0

TL = TLMat(rand(n,r), rand(n,r) + 1i * rand(n,r));
sigma = rand + 1i * rand;
TL_true = full(TL) + diag(sigma * ones(n,1));
TL = TL + sigma * toepeye(n);
nTL = norm(full(TL));
testCase.assertEqual(class(TL), 'TLMat');
testCase.assertEqual(full(TL), TL_true, 'AbsTol', 2*nTL*eps, 'RelTol', 2*nTL*eps);
testCase.assertEqual(drank(TL), r + 1); % Fails with prob. 0

% Sum of random matrices
n = 11;
dr = 3;
TL = TLMat(rand(n,dr), rand(n,dr) + 1i * rand(n,dr));
[c,r, T] = random_toeplitz(n,n);
TM = ToepMat(c,r);

nfact = 4 * (norm(full(TM)) + norm(full(TL)));
testCase.assertEqual(class(TL + TM), 'TLMat');
testCase.assertEqual(class(TL - TM), 'TLMat');
testCase.assertEqual(class(TM + TL), 'TLMat');
testCase.assertEqual(class(TM - TL), 'TLMat');
testCase.assertEqual(full(TL + TM), full(TL) + T, ...
    'RelTol', nfact*eps, 'AbsTol', nfact*eps);
testCase.assertEqual(full(TL - TM), full(TL) - T, ...
    'RelTol', nfact*eps, 'AbsTol', nfact*eps);
testCase.assertEqual(full(TM + TL), T + full(TL), ...
    'RelTol', nfact*eps, 'AbsTol', nfact*eps);
testCase.assertEqual(full(TM - TL), T - full(TL), ...
    'RelTol', nfact*eps, 'AbsTol', nfact*eps);

% These four tests fail with probability 0
testCase.assertEqual(drank(TL + TM), dr+2);
testCase.assertEqual(drank(TL - TM), dr+2);
testCase.assertEqual(drank(TM + TL), dr+2);
testCase.assertEqual(drank(TM - TL), dr+2);

end

function test_unary_minus(testCase)

TL = TLMat([]);
testCase.assertTrue(isempty(full(TL - TL)));

TL = tleye(1);
testCase.assertEqual(full(-TL), -1);

TL = TLMat(randn(9,3), randn(9,3) + 1i*randn(9,3));
testCase.assertEqual(full(-TL), -full(TL));
end

function test_transpose(testCase)

TL = TLMat([]);
testCase.assertEqual(class(TL.'), 'TLMat');
testCase.assertTrue(isempty(full(TL.')));

TL = TLMat(-1i);
testCase.assertEqual(class(TL.'), 'TLMat');
testCase.assertEqual(full(TL.'), full(TL).');
testCase.assertEqual(full((TL.').'), full(TL));

TL = tleye(6);
testCase.assertEqual(class(TL.'), 'TLMat');
testCase.assertEqual(full(TL.'), full(TL).', 'AbsTol', 4*eps, 'RelTol', 4*eps);
testCase.assertEqual(full((TL.').'), full(TL), 'AbsTol', 4*eps, 'RelTol', 4*eps);

TL = TLMat(randn(9,4), 1i * rand(9,4));
nfact = 16*norm(full(TL));
testCase.assertEqual(class(TL.'), 'TLMat');
testCase.assertEqual(full(TL.'), full(TL).', ...
    'Abstol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(full((TL.').'), full(TL), ...
    'Abstol', nfact*eps, 'RelTol', nfact*eps);

end

function test_ctranspose(testCase)

TL = TLMat([]);
testCase.assertEqual(class(TL'), 'TLMat');
testCase.assertTrue(isempty(full(TL')));

TL = TLMat(-1i);
testCase.assertEqual(class(TL'), 'TLMat');
testCase.assertEqual(full(TL'), full(TL)');
testCase.assertEqual(full((TL')'), full(TL));

TL = tleye(6);
testCase.assertEqual(class(TL'), 'TLMat');
testCase.assertEqual(full(TL'), full(TL)', 'AbsTol', 4*eps, 'RelTol', 4*eps);
testCase.assertEqual(full((TL')'), full(TL), 'AbsTol', 4*eps, 'RelTol', 4*eps);

TL = TLMat(randn(9,4), 1i * rand(9,4));
nfact = 16*norm(full(TL));
testCase.assertEqual(class(TL'), 'TLMat');
testCase.assertEqual(full(TL'), full(TL)', ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(full((TL')'), full(TL), ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
end

function test_transpose_rank(testCase)
% A curiosity of the Zp/Zm Sylvester operator is that the rank may increase
% under transposition.   For TL matrices that are in fact Toeplitz
% matrices, or "derive" from such a matrix, the rank should not increase.
n = 13;
[c,r] = random_toeplitz(n,n);
[G, B] = toepgen(c,r);
T = TLMat(G, B);

testCase.assertEqual(drank(T), 2);
testCase.assertEqual(drank(T'), 2); % This will fail sometimes with drank=3

G = randn(n,4);
B = randn(n,4);
T = TLMat(G,B);
testCase.assertEqual(drank(T), 4);
testCase.assertEqual(drank(T'), 6);
end

function test_mtimes_scalar(testCase)
TL = TLMat([]);
s = 8;
testCase.assertTrue(isempty(full(s * TL)));

TL = TLMat(1i);
s = 1i;
testCase.assertEqual(class(s * TL), 'TLMat');
testCase.assertEqual(class(TL * s), 'TLMat');
testCase.assertEqual(full(s * TL), -1);
testCase.assertEqual(full(TL * s), -1);

TL = tleye(5);
s = -2;
testCase.assertEqual(class(s * TL), 'TLMat');
testCase.assertEqual(class(TL * s), 'TLMat');
testCase.assertEqual(full(s * TL), -2 * eye(5));
testCase.assertEqual(full(TL * s), -2 * eye(5));

[c,r,T] = random_toeplitz(7,7);
TL = TLMat(c,r);
s = 1i * pi;
testCase.assertEqual(class(s * TL), 'TLMat');
testCase.assertEqual(class(TL * s), 'TLMat');
testCase.assertEqual(full(s * TL), s * T, 'AbsTol', 30*eps, 'RelTol', 30*eps);
testCase.assertEqual(full(TL * s), s * T, 'AbsTol', 30*eps, 'RelTol', 30*eps);

TL = TLMat(1i * rand(12,4), rand(12,4));
nTL = norm(full(TL));
s = -2 + 4i;
testCase.assertEqual(class(s * TL), 'TLMat');
testCase.assertEqual(class(TL * s), 'TLMat');
testCase.assertEqual(full(s * TL), s * full(TL), ...
    'AbsTol', 2*nTL*eps, 'RelTol', 2*nTL*eps);
testCase.assertEqual(full(TL * s), s * full(TL), ...
    'AbsTol', 2*nTL*eps, 'RelTol', 2*nTL*eps);

end

function test_mtimes_toepmat(testCase)

TL = tleye(5);
TM = toepeye(4);
testCase.assertError( @() TL * TM, 'tlcomp:InconsistentInput');
testCase.assertError( @() TM * TL, 'tlcomp:InconsistentInput');

TL = TLMat([]);
TM = ToepMat([], []);
testCase.assertEqual(class(TL * TM), 'TLMat');
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertTrue(isempty(full(TL*TM)));
testCase.assertTrue(isempty(full(TM*TL)));

TL = TLMat(1i);
TM = ToepMat(1i,1i);
testCase.assertEqual(class(TL * TM), 'TLMat');
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertEqual(full(TL * TM), -1);
testCase.assertEqual(full(TM * TL), -1);

TL = tleye(6);
TM = toepeye(6);
testCase.assertEqual(class(TL * TM), 'TLMat');
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertEqual(full(TL * TM), eye(6), 'AbsTol', eps);
testCase.assertEqual(full(TM * TL), eye(6), 'Abstol', eps);
testCase.assertEqual(drank(TM * TL), 1);

TL = TLMat(rand(12,3), rand(12,3) + 1i * rand(12,3));
[c,r,T] = random_toeplitz(12,12);
TM = ToepMat(c,r);
nfact = 4 * (norm(full(TM)) * norm(full(TL)));
testCase.assertEqual(class(TL * TM), 'TLMat');
testCase.assertEqual(class(TM * TL), 'TLMat');
testCase.assertEqual(full(TL * TM), full(TL) * T, ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(full(TM * TL), T * full(TL), ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(drank(TM * TL), 5);
end

function test_mtimes_tlmat(testCase)

TL1 = tleye(8);
TL2 = TLMat([1,2,3], [1, -2, -1]);
testCase.assertError( @() TL1 * TL2, 'tlcomp:InconsistentInput');
testCase.assertError( @() TL2 * TL1, 'tlcomp:InconsistentInput');

TL1 = TLMat(0); % 1-by-1 drank 0
TL2 = TLMat(0);
testCase.assertEqual(size(TL1 * TL2), [1,1]);
testCase.assertEqual(full(TL1 * TL2), 0);

TL1 = tleye(1);
TL2 = TLMat(3i);
testCase.assertEqual(class(TL1 * TL2), 'TLMat');
testCase.assertEqual(class(TL2 * TL1), 'TLMat');
testCase.assertEqual(full(TL1 * TL2), 3i);
testCase.assertEqual(full(TL2 * TL1), 3i);
testCase.assertEqual(drank(TL1 * TL2), 1);

TL1 = TLMat(rand(8,2), 1i*rand(8,2));
TL2 = tleye(8);
testCase.assertEqual(class(TL1 * TL2), 'TLMat');
testCase.assertEqual(class(TL2 * TL1), 'TLMat');
testCase.assertEqual(full(TL1 * TL2), full(TL1), ...
    'AbsTol', 100*eps, 'RelTol', 100*eps);
testCase.assertEqual(full(TL2 * TL1), full(TL1), ...
    'AbsTol', 50*eps, 'RelTol', 50*eps);

TL1 = TLMat(rand(8,3), 1i*rand(8,3));
TL2 = TLMat(ones(8,1), ones(8,1));
n1 = norm(full(TL1));
n2 = norm(full(TL2));
testCase.assertEqual(class(TL1 * TL2), 'TLMat');
testCase.assertEqual(class(TL2 * TL1), 'TLMat');
testCase.assertEqual(full(TL1 * TL2), full(TL1) * ones(8), ...
    'AbsTol', 2*n1*n2*eps, 'RelTol', 2*n1*n2*eps);
testCase.assertEqual(full(TL2 * TL1), ones(8) * full(TL1), ...
    'AbsTol', n1*n2*eps, 'RelTol', n1*n2*eps);

TL1 = TLMat(rand(9,1), 1i*rand(9,1), 'GB');
TL2 = TLMat(1i*randn(9,1), randn(9,1), 'GB');
testCase.assertEqual(class(TL1 * TL2), 'TLMat');
testCase.assertEqual(class(TL2 * TL1), 'TLMat');
testCase.assertEqual(full(TL1 * TL2), full(TL1) * full(TL2), ...
    'AbsTol', 500*eps, 'RelTol', 500*eps);
testCase.assertEqual(full(TL2 * TL1), full(TL2) * full(TL1), ...
    'AbsTol', 500*eps, 'RelTol', 500*eps);
testCase.assertEqual(drank(TL1 * TL2), 3);



TL1 = TLMat(rand(9,3), rand(9,3));
TL2 = TLMat(1i*randn(9,3), randn(9,3));
testCase.assertEqual(class(TL1 * TL2), 'TLMat');
testCase.assertEqual(class(TL2 * TL1), 'TLMat');
testCase.assertEqual(full(TL1 * TL2), full(TL1) * full(TL2), ...
    'AbsTol', 1e4*eps, 'RelTol', 1e4*eps);
testCase.assertEqual(full(TL2 * TL1), full(TL2) * full(TL1), ...
    'AbsTol', 1e4*eps, 'RelTol', 1e4*eps);

testCase.assertEqual(drank(TL1 * TL2), 7);

end

function test_mtimes_double_matrix(testCase)

TL = tleye(8);
A = rand(7);
testCase.assertError( @() TL * A, 'tlcomp:InconsistentInput');
testCase.assertError( @() A * TL, 'tlcomp:InconsistentInput');

TL = TLMat([]);
A = [];
testCase.assertEqual(class(TL * A), 'double');
testCase.assertEqual(class(A * TL), 'double');
testCase.assertEqual(TL * A, []);
testCase.assertEqual(A * TL, []);

TL = tleye(3);
A = eye(3);
testCase.assertEqual(class(TL * A), 'double');
testCase.assertEqual(class(A * TL), 'double');
testCase.assertEqual(TL * A, eye(3), 'AbsTol', 2*eps);
testCase.assertEqual(A * TL, eye(3), 'AbsTol', 2*eps);

TL = TLMat(rand(3), rand(3));
nTL = norm(full(TL));
A = rand(3);
nA = norm(A);

nfact = 16 * nTL * nA;

testCase.assertEqual(class(TL * A), 'double');
testCase.assertEqual(class(A * TL), 'double');
testCase.assertEqual(TL * A, full(TL) * A, ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(A * TL, A * full(TL), ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);

TL = TLMat(randn(12,3) + 1i * randn(12,3), randn(12,3));
A = 1i * randn(12);
nA = norm(A);
nTL = norm(full(TL));

nfact = 4 * nTL * nA;

testCase.assertEqual(class(TL * A), 'double');
testCase.assertEqual(class(A * TL), 'double');
testCase.assertEqual(TL * A, full(TL) * A, ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
testCase.assertEqual(A * TL, A * full(TL), ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
end

function test_mtimes_double_vector(testCase)

TL = tleye(9);
x = ones(9,1);
testCase.assertEqual(TL * x, x, 'AbsTol', eps);

TL = TLMat(randn(9,3), randn(9,3));
nTL = norm(full(TL));
x = randn(9,1);
testCase.assertEqual(TL * x, full(TL) * x, 'AbsTol', 4*nTL*eps, 'RelTol', 4*nTL*eps);
end

function test_mldivide_double(testCase)

TL = TLMat(randn(14,4), randn(14,4));
b = randn(14,3);
x = TL \ b;
xx = full(TL) \ b;

nfact = cond(full(TL)) * norm(b);
testCase.assertTrue(isreal(x));
testCase.assertEqual(x, xx, 'AbsTol', nfact*eps, 'RelTol', nfact*eps);

end

function test_mldivide_tlmat(testCase)

TL1 = tleye(5);
TL2 = tleye(5);
testCase.assertEqual(class(TL1 \ TL2), 'TLMat');
testCase.assertEqual(drank(TL1 \ TL2), 1);
testCase.assertEqual(full(TL1 \ TL2), eye(5), 'AbsTol', eps);

TL1 = TLMat(randn(12,3), 1i * randn(12,3)) + tleye(12);
TL2 = TLMat(1i * randn(12,2), randn(12,2));
D = TL1 \ TL2;
nfact = cond(full(TL1)) * norm(full(TL2));
testCase.assertEqual(drank(D), 6);
testCase.assertEqual(full(D), full(TL1) \ full(TL2), ...
    'AbsTol', nfact*eps, 'RelTol', nfact*eps);
end

function test_mldivide_toepmat(testCase)

TL = tleye(9);
[c,r,T] = random_toeplitz(9,9);
TM = ToepMat(c,r);
D = TL\TM;
testCase.assertTrue(isa(D, 'TLMat'));
testCase.assertEqual(full(D), T, 'AbsTol', 100*eps, 'RelTol', 100*eps);
testCase.assertEqual(drank(D), 2);

TL = TLMat(randn(9,4) + 1i*randn(9,4), randn(9,4) + 1i * randn(9,4));
D = TL \ TM;
testCase.assertEqual(full(D), full(TL) \ T, ...
    'AbsTol', 1e4*eps, 'RelTol', 1e4*eps);
testCase.assertEqual(drank(D), 7);

end

function test_mrdivide_scalar(testCase)
T = TLMat(1,1);
s = 1;
T = T/s;
testCase.assertEqual(full(T), 1);
testCase.assertTrue(isa(T, 'TLMat'));

s = 2;
T = T/s;
testCase.assertEqual(full(T), .5);

s = .5i;
T = T/s;
testCase.assertEqual(full(T), -1i);
end


function test_mrdivide_toepmat(testCase)
TL = TLMat(2i,2i);
B = ToepMat(3,3);
testCase.assertTrue(isa(B/TL, 'TLMat'));
testCase.assertEqual(full(B/TL), -3i/2);

TL = tleye(9);
B = ToepMat(1:9);
testCase.assertTrue(isa(B/TL, 'TLMat'));
testCase.assertEqual(full(B/TL), full(B), 'AbsTol', 64*eps, 'RelTol', 64*eps);

TL = TLMat(randn(7,2), randn(7,2));
B = ToepMat(randn(7,1));
nfact = 16 * cond(full(TL)) * norm(full(B));
testCase.assertEqual(full(B/TL), ...
    full(B)/full(TL), 'AbsTol', nfact*eps, 'RelTol', nfact*eps);
end


function test_mrdivide_double(testCase)

TL = tleye(9);
B = reshape(1:27, 3, 9);
testCase.assertTrue(isa(B/TL, 'double'));
testCase.assertEqual(B/TL, B, 'AbsTol', 32*eps, 'RelTol', 32*eps);

B = eye(9);
testCase.assertTrue(isa(B/TL, 'double'));
testCase.assertEqual(B/TL, eye(9), 'AbsTol', 4*eps, 'RelTol', 4*eps);

TL = TLMat(randn(9,4), randn(9,4));
b = randn(1,9) + 1i * randn(1,9);
nfact = cond(full(TL)) * norm(b);
testCase.assertEqual(b/TL, b/full(TL), 'AbsTol', nfact*eps, 'RelTol', nfact*eps);

testCase.assertError(@() b' / TL, 'tlcomp:InconsistentInput');
end

function test_mrdivide_tlmat(testCase)
TL = TLMat(2i,2i);
B = TLMat(3,3);
testCase.assertTrue(isa(B/TL, 'TLMat'));
testCase.assertEqual(full(B/TL), -3i/2);

TL = tleye(8);
B = TLMat(randn(8,2), randn(8,2));
nfact = 16*norm(full(B));
testCase.assertTrue(isa(B/TL, 'TLMat'));
testCase.assertEqual(full(B/TL), full(B)/full(TL), ...
    'Abstol', nfact*eps, 'RelTol', nfact*eps);

TL = TLMat(rand(12,3), 1i * randn(12,3)) + tleye(12);
B = TLMat(randn(12,3), randn(12,3));
nfact = 16 * cond(full(TL)) * norm(full(B));
testCase.assertEqual(full(B/TL), full(B)/full(TL), ...
    'RelTol', nfact * eps, 'AbsTol', nfact*eps);
end

function test_norm_eye(testCase)
n = 14;
TL = tleye(14);
testCase.assertEqual(norm(TL, 1), 1);
testCase.assertEqual(norm(TL, inf), 1);
testCase.assertEqual(norm(TL, 'inf'), 1);
testCase.assertEqual(norm(TL, 'fro'), sqrt(n));
end

function test_norm_toeplitz(testCase)

[c, r, T] = random_toeplitz(7,7);
TL = TLMat(c,r);
testCase.assertEqual(norm(TL, 1), norm(T, 1), 'RelTol', 4*eps, 'AbsTol', 4*eps);
testCase.assertEqual(norm(TL, 'inf'), norm(T, 'inf'), 'RelTol', 4*eps, 'AbsTol', 4*eps);
testCase.assertEqual(norm(TL, inf), norm(T, inf), 'RelTol', 4*eps, 'AbsTol', 4*eps);
testCase.assertEqual(norm(TL, 'fro'), norm(T, 'fro'), 'RelTol', 4*eps, 'AbsTol', 4*eps);

testCase.assertWarning(@() norm(TL), 'tlcomp:Unsupported');
testCase.assertWarning(@() norm(TL, 2), 'tlcomp:Unsupported');
testCase.assertError(@() norm(TL, 3), 'tlcomp:InconsistentInput');

warning('OFF', 'tlcomp:Unsupported');
testCase.assertEqual(norm(TL), norm(T), 'RelTol', 8*eps, 'AbsTol', 8*eps);
testCase.assertEqual(norm(TL, 2), norm(T, 2), 'RelTol', 8*eps, 'AbsTol', 8*eps);
warning('ON', 'tlcomp:Unsupported');
end

function test_norm_randomtl(testCase)
TL = TLMat(randn(6,3), randn(6,3) + 1i * rand(6,3));
T = full(TL);
testCase.assertEqual(norm(TL, 1), norm(T, 1), 'RelTol', 4*eps, 'AbsTol', 4*eps);
testCase.assertEqual(norm(TL, 'inf'), norm(T, 'inf'), 'RelTol', 4*eps, 'AbsTol', 4*eps);
testCase.assertEqual(norm(TL, inf), norm(T, inf), 'RelTol', 4*eps, 'AbsTol', 4*eps);
testCase.assertEqual(norm(TL, 'fro'), norm(T, 'fro'), 'RelTol', 4*eps, 'AbsTol', 4*eps);

warning('OFF', 'tlcomp:Unsupported');
testCase.assertEqual(norm(TL), norm(T));
testCase.assertEqual(norm(TL, 2), norm(T, 2));
warning('ON', 'tlcomp:Unsupported');
end

function test_mpower_monomial(testCase)
TM = TLMat(0);
testCase.assertEqual(full(TM^0), 1);
testCase.assertTrue(isa(TM^0, 'TLMat'));
for s = 1:3
    testCase.assertEqual(full(TM^s), 0);
    testCase.assertTrue(isa(TM^s, 'TLMat'));
end

TM = tleye(1);
for s = 0:4
    testCase.assertEqual(full(TM^s), 1.0);
end

TM = tleye(3);
for s = 0:4
    testCase.assertEqual(full(TM^s), eye(3), 'AbsTol', eps);
    testCase.assertTrue(isa(TM^s, 'TLMat'));

end

TM = TLMat(1i*randn(16,3), randn(16,3));
T = full(TM);
for s = 0:4
    Tpow_true = T^s;
    Tpow = TM^s;
    nT = norm(Tpow_true);
    testCase.assertEqual(full(Tpow), Tpow_true, ...
        'AbsTol', 4*nT*eps, 'RelTol', 4*nT*eps);
end

end

function test_gennorm(testCase)
TM = TLMat(0);
testCase.assertEqual(TM.gennorm(), 0.0);

TM = tleye(6);
testCase.assertEqual(TM.gennorm(), 2.0);

TM = TLMat(randn(7,3), randn(7,3));
testCase.assertEqual(TM.gennorm(), norm(TM.G * TM.B'), 'RelTol', 4*eps);
end

