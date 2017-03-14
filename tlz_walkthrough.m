%% tlzstein -- Toeplitz and Toeplitz-like matrices w.r.t. Z-Stein displacement
%
% This collection of functions and classes contains a number of tools to
% carry out matrix computations with Toeplitz and Toeplitz-like matrices
% with respect to the Stein displacement operator
%
%   D(T) := T - Z * T * Z'
%
% where Z is the _downshift_ matrix, viz. zeros(n) + diag(ones(n-1,1),-1).
% Any matrix T such that D(T) has _low rank_ is said to be Toeplitz-like
% (TL).

%% Setup path
addpath('matlab')  % Code is sitting here
addpath('drsolve') % Used for GKO


%% Basic operations
%
% Run these examples youself!

% Generate data for two toeplitz matrices
[c1, r1] = random_toeplitz(8, 8);
[c2, r2] = random_toeplitz(8, 8);

% We provide a class |ToepMat|
TM1 = ToepMat(c1, r1);
TM2 = ToepMat(c2, r2);

% Addition, scalar multiplication yield a |ToepMat|
disp(TM1 + TM2)
disp(TM1 - TM2)
disp(2i*pi * TM1)

% Mixed |ToepMat| / double array products
disp(TM1 * ones(8,1));
disp(TM1 * ones(8,8));

% Product of two Toeplitz matrices yields a |TLMat| object, d-rank is 4
disp(TM1 * TM2)

% Convert to full double matrix
disp(full(TM1))

% Evaluate deg-6 Taylor polynomial of expm(TM1)
p = 1./factorial(6:-1:0);
E = polyvalm(p, TM1);
disp(E); % Result is a TLMat
EE = polyvalm(p, full(TM1)); % Compare with dense computation
disp(norm(E - EE, 'fro') / norm(EE, 'fro'));

% Toeplitz-inverse is a |TLMat| of d-rank 2
disp(inv(TM1));

% Shifted inverseXS, useful for partial fraction expansion
sigma = exp(3i/4 * pi);
disp(inv(TM1 - sigma * toepeye(8)));

% Solve a linear system
b = TM1 * ones(8,1);
x = TM1 \ b;
disp(norm(x - ones(8,1)));

% It's also quite fast, using the GKO algorithm
[c,r,T] = random_toeplitz(4096,4096);
TMbig = ToepMat(c,r);
b = randn(4096,1);
tic; x = TMbig \ b; toc
tic; xx = T\b; toc
disp(norm(x - xx)/norm(xx));

% Rational function evaluation, here: Cayley transform
e2 = zeros(16,1);
e2(2) = 1;
TM = ToepMat(-e2, e2); % Skew-symm
rT = TM.ratevalm([1,1], [1,-1]); % This is a |TLMat|
disp(norm(full(rT) - ( (full(TM) - eye(16)) \ (full(TM) + eye(16) ) )));

% Can also be computed with \
rT = (TM - toepeye(16)) \ (TM + toepeye(16)); % Again a |TLMat|
disp(norm(full(rT) - ( (full(TM) - eye(16)) \ (full(TM) + eye(16) ) )));





