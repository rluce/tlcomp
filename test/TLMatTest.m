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
A = 6;
TL = TLMat(A);
testCase.assertEqual(TL.G * TL.B', A);

G = 2;
B = 3;
TL = TLMat(G,B, 'GB');
testCase.assertEqual(TL.G * TL.B', A);

G = [1, 1];
B = [3, 3];
TL = TLMat(G,B, 'GB');
testCase.assertEqual(TL.G * TL.B', A);

c = 6;
r = 6;
TL = TLMat(c,r);
testCase.assertEqual(TL.G * TL.B', A);

c = 6;
r = 1;
testCase.assertError( @() TLMat(c,r), 'tlzstein:InconsistentInput');


end


function test_construct_eye(testCase)
n = 8;

E = eye(8);
e1 = E(:,1);
z = zeros(n,1);

TL = TLMat(E);
testCase.assertEqual(TL.G * TL.B', e1 * e1');

TL = TLMat(e1,e1);
testCase.assertEqual(TL.G * TL.B', e1 * e1');

TL = TLMat(e1,e1, 'GB');
testCase.assertEqual(TL.G * TL.B', e1 * e1');

TL = TLMat([e1,z], [e1,z], 'GB');
testCase.assertEqual(TL.G * TL.B', e1 * e1');

TL = TLMat([e1,z], [e1,z]);
testCase.assertEqual(TL.G * TL.B', e1 * e1');

end

function test_construct_random(testCase)
n = 9;
[c,r,T] = random_toeplitz(n,n);

TL = TLMat(T);
T2 = stein_reconstruction(TL.G, TL.B);
testCase.assertEqual(T2, T, 'RelTol', 100*eps);

TL = TLMat(c,r);
T2 = stein_reconstruction(TL.G, TL.B);
testCase.assertEqual(T2, T, 'RelTol', 100*eps);

[G, B] = stein_generator(c,r);
TL = TLMat(G,B);
T2 = stein_reconstruction(TL.G, TL.B);
testCase.assertEqual(T2, T, 'RelTol', 100*eps);

GG = [2 * G, -G];
BB = [B, B];
TL = TLMat(GG, BB);
T2 = stein_reconstruction(TL.G, TL.B);
testCase.assertEqual(T2, T, 'RelTol', 100*eps);


end

function test_construct_badinput(testCase)

testCase.assertError( @() TLMat(rand(4,2), rand(4,3)), ...
    'tlzstein:InconsistentInput');

end


function test_size(testCase)

TL = TLMat([]);
testCase.assertEqual(size(TL), [0,0]);

TL = TLMat(1);
testCase.assertEqual(size(TL), [1,1]);

TL = TLMat([1,2,5], [1,2,5]);
testCase.assertEqual(size(TL), [3,3]);

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
testCase.assertEqual(full(TLp1), full(TL) + 1, 'RelTol', 100*eps);

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

% Could be represented as TL, but for consistency we want dense
B = A + TL;
testCase.assertEqual(class(B), 'double');
testCase.assertTrue(isempty(B));

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

B = TL + zeros(n);
testCase.assertEqual(class(B), 'TLMat');
testCase.assertEqual(drank(B), 3);
testCase.assertEqual(full(B), full(TL), 'RelTol', 10*eps);

sigma = rand(1,1) + 1i * rand(1,1);
B = TL + sigma * eye(n);
testCase.assertEqual(class(B), 'TLMat');
testCase.assertEqual(drank(B), 4); % fails with prob 0
testCase.assertEqual(full(B), full(TL) + sigma * eye(n), 'RelTol', 10*eps);

[~, ~, T] = random_toeplitz(n,n);
testCase.assertEqual(class(TL + T), 'TLMat');
testCase.assertEqual(class(TL - T), 'TLMat');
testCase.assertEqual(class(T + TL), 'TLMat');
testCase.assertEqual(class(T - TL), 'TLMat');

testCase.assertEqual(drank(TL + T), 5); % fails with prob 0
testCase.assertEqual(drank(TL - T), 5); % fails with prob 0
testCase.assertEqual(drank(T + TL), 5); % fails with prob 0
testCase.assertEqual(drank(T - TL), 5); % fails with prob 0

testCase.assertEqual(full(TL + T), full(TL) + T, 'RelTol', 100*eps);
testCase.assertEqual(full(TL - T), full(TL) - T, 'RelTol', 100*eps);
testCase.assertEqual(full(T + TL), T + full(TL), 'RelTol', 100*eps);
testCase.assertEqual(full(T - TL), T - full(TL), 'RelTol', 100*eps);

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
    'RelTol', 100*eps);

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
testCase.assertEqual(class(TL), 'TLMat');
testCase.assertEqual(full(TL), TL_true, 'RelTol', 10*eps);
testCase.assertEqual(drank(TL), r + 1); % Fails with prob. 0

TL = TLMat(rand(n,r), rand(n,r) + 1i * rand(n,r));
sigma = rand + 1i * rand;
TL_true = full(TL) + diag(sigma * ones(n,1));
TL = TL + sigma * toepeye(n);
testCase.assertEqual(class(TL), 'TLMat');
testCase.assertEqual(full(TL), TL_true);
testCase.assertEqual(drank(TL), r + 1); % Fails with prob. 0

% Sum of random matrices
n = 11;
r = 3;
TL = TLMat(rand(n,r), rand(n,r) + 1i * rand(n,r));
[c,r, T] = random_toeplitz(n,n);
TM = ToepMat(c,r);
testCase.assertEqual(class(TL + TM), 'TLMat');
testCase.assertEqual(class(TL - TM), 'TLMat');
testCase.assertEqual(class(TM + TL), 'TLMat');
testCase.assertEqual(class(TM - TL), 'TLMat');
testCase.assertEqual(full(TL + TM), full(TL) + T);
testCase.assertEqual(full(TL - TM), full(TL) - T);
testCase.assertEqual(full(TM + TL), T + full(TL));
testCase.assertEqual(full(TM - TL), T - full(TL));

% These four tests fail with probability 0
testCase.assertEqual(drank(TL + TM), r+2);
testCase.assertEqual(drank(TL - TM), r+2);
testCase.assertEqual(drank(TM + TL), r+2);
testCase.assertEqual(drank(TM - TL), r+2);

end

function test_unary_minus(testCase)

TL = TLMat([]);
testCase.assertTrue(isempty(full(TL - TL)));

TL = tleye(1);
testCase.assertEqual(full(-TL), -1);

TL = TLMat(randn(9,3), randn(9,3) + 1i*randn(9,3));
testCase.assertEqual(full(-TL), -full(TL));
end

