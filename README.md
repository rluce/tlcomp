TLComp -- Computing with Toeplitz-like matrices in Matlab
======

Toeplitz matrices, and matrices arising from arithmetic operations
among Toeplitz matrices enjoy a _low rank representation_ w.r.t. to
so-called displacement operators [1].  This toolbox offers
convenient operations with such strucutred matices while taking
advantage of such a low rank representation.

Installation
---

1. TLComp relies on the "drsolve" [2] package for operations that involve
solving systems of linear equations, please download and install it from http://bugs.unica.it/~gppe/software/drsolve/drsolve1.0.zip

2. Add the 'matlab' folder below this directory to the matlab path,
   viz. `addpath(fullfile(pwd, 'matlab'))`

Examples
---

```
%% Generate data for two toeplitz matrices
[c1, r1] = random_toeplitz(8, 8);
[c2, r2] = random_toeplitz(8, 8);

%% We provide a class |ToepMat|
TM1 = ToepMat(c1, r1);
TM2 = ToepMat(c2, r2);

%% Addition, scalar multiplication yield a |ToepMat|
disp(TM1 + TM2)
disp(TM1 - TM2)
disp(2i*pi * TM1)

%% Mixed |ToepMat| / double array products
disp(TM1 * ones(8,1));
disp(TM1 * ones(8,8));

%% Product of two Toeplitz matrices yields a |TLMat| object, d-rank is 4
disp(TM1 * TM2)

%% Convert to full double matrix
disp(full(TM1))

%% Evaluate deg-6 Taylor polynomial of expm(TM1)
p = 1./factorial(6:-1:0);
E = polyvalm(p, TM1);
disp(E); % Result is a TLMat
EE = polyvalm(p, full(TM1)); % Compare with dense computation
disp(norm(E - EE, 'fro') / norm(EE, 'fro'));

%% Toeplitz-inverse is a |TLMat| of d-rank 2
disp(inv(TM1));

%% Shifted inverse, useful for partial fraction expansion
sigma = exp(3i/4 * pi);
disp(inv(TM1 - sigma * toepeye(8)));

%% Solve a linear system
b = TM1 * ones(8,1);
x = TM1 \ b;
disp(norm(x - ones(8,1)));

%% It's also quite fast, using the GKO algorithm
[c,r,T] = random_toeplitz(4096,4096);
TMbig = ToepMat(c,r);
b = randn(4096,1);
tic; x = TMbig \ b; toc
tic; xx = T\b; toc
disp(norm(x - xx)/norm(xx));

%% Rational function evaluation, here: Cayley transform
e2 = zeros(16,1);
e2(2) = 1;
TM = ToepMat(-e2, e2); % Skew-symm
rT = TM.ratevalm([1,1], [1,-1]); % This is a |TLMat|
disp(norm(full(rT) - ( (full(TM) - eye(16)) \ (full(TM) + eye(16) ) )));

%% Can also be computed with \
rT = (TM - toepeye(16)) \ (TM + toepeye(16)); % Again a |TLMat|
disp(norm(full(rT) - ( (full(TM) - eye(16)) \ (full(TM) + eye(16) ) )));

%% Determinant is supported through GKO/LU factorization
TM = ToepMat([1,2,3,4]);
disp(det(TM))

%% Construct A TLMat object from a generator
G = orth(randn(8,4));
B = orth(randn(8,4)) * diag([10,1,0.01, 1e-12]);
TL1 = TLMat(G, B);  % drank is 4
TL2 = TL1.truncate_tol(1e-8); % remove ranks relatively smaller than 1e-8
fprintf('Trunctated from %d to %d, generator difference: %.2e\n', ...
    drank(TL1), drank(TL2), gennorm(TL1 - TL2));
```

More background
---

The displacement operator we use in this toolbox is the Sylvester-type operator
```
   D(A) = Z_{1} A - A * Z_{-1}
```
where `Z_s` is the s-circulant matrix, i.e., for n=4
```
Z_s = [
    0 0 0 s
    1 0 0 0
    0 1 0 0
    0 0 1 0
]
```



References
----


[1] Kailath, T., & Sayed, A. H. (Eds.). (1999). Fast reliable algorithms for matrices with structure. Society for Industrial and Applied Mathematics (SIAM), Philadelphia, PA. https://doi.org/10.1137/1.9781611971354

[2] A. Aricò and G. Rodriguez.  A fast solver for linear systems with
displacement structure.  Numer. Algorithms, 55(4):529-556, 2010.  DOI:
https://doi.org/10.1007/s11075-010-9421-x
