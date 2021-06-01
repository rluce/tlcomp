TLComp -- Computing with Toeplitz-like matrices in Matlab
======

Toeplitz matrices, and matrices arising from arithmetic operations
among Toeplitz matrices enjoy a _low rank representation_ w.r.t. to
so-called displacement operators [1].  This toolbox offers
convenient and fast operations with such structured matices, taking
advantage of their a low rank representation.

Installation
---

1. TLComp relies on the "drsolve" [2] package for operations that involve
solving systems of linear equations, please download and install it from http://bugs.unica.it/~gppe/software/drsolve/drsolve1.0.zip

3. Download and unpack the [TLComp toolbox](https://github.com/rluce/tlcomp/releases/download/v0.1.0/tlcomp-v0.1.0.tar.gz)

2. Within Matlab add the installation folder to the matlab path,
   viz. `addpath(fullfile(pwd, 'tlcomp-v0.1.0'))`.

Examples
---

```matlab
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
```
```
  8x8 ToepMat
  8x8 ToepMat
  8x8 ToepMat
```
```matlab
%% Mixed |ToepMat| / double array products yield full arrays
disp(TM1 * ones(8,1));
```
```
  -0.8618 - 0.8925i
   0.0425 - 0.7769i
  -1.4210 - 0.6578i
  -0.0681 - 0.5681i
   1.4395 - 1.3148i
  -0.2738 - 2.0862i
   0.8891 - 1.8790i
  -2.9664 - 2.8394i
```
```matlab
disp(TM1 * ones(8,8));
```
```
  Columns 1 through 5

  -0.8618 - 0.8925i  -0.8618 - 0.8925i  -0.8618 - 0.8925i  -0.8618 - 0.8925i  -0.8618 - 0.8925i
   0.0425 - 0.7769i   0.0425 - 0.7769i   0.0425 - 0.7769i   0.0425 - 0.7769i   0.0425 - 0.7769i
  -1.4210 - 0.6578i  -1.4210 - 0.6578i  -1.4210 - 0.6578i  -1.4210 - 0.6578i  -1.4210 - 0.6578i
  -0.0681 - 0.5681i  -0.0681 - 0.5681i  -0.0681 - 0.5681i  -0.0681 - 0.5681i  -0.0681 - 0.5681i
   1.4395 - 1.3148i   1.4395 - 1.3148i   1.4395 - 1.3148i   1.4395 - 1.3148i   1.4395 - 1.3148i
  -0.2738 - 2.0862i  -0.2738 - 2.0862i  -0.2738 - 2.0862i  -0.2738 - 2.0862i  -0.2738 - 2.0862i
   0.8891 - 1.8790i   0.8891 - 1.8790i   0.8891 - 1.8790i   0.8891 - 1.8790i   0.8891 - 1.8790i
  -2.9664 - 2.8394i  -2.9664 - 2.8394i  -2.9664 - 2.8394i  -2.9664 - 2.8394i  -2.9664 - 2.8394i

  Columns 6 through 8

  -0.8618 - 0.8925i  -0.8618 - 0.8925i  -0.8618 - 0.8925i
   0.0425 - 0.7769i   0.0425 - 0.7769i   0.0425 - 0.7769i
  -1.4210 - 0.6578i  -1.4210 - 0.6578i  -1.4210 - 0.6578i
  -0.0681 - 0.5681i  -0.0681 - 0.5681i  -0.0681 - 0.5681i
   1.4395 - 1.3148i   1.4395 - 1.3148i   1.4395 - 1.3148i
  -0.2738 - 2.0862i  -0.2738 - 2.0862i  -0.2738 - 2.0862i
   0.8891 - 1.8790i   0.8891 - 1.8790i   0.8891 - 1.8790i
  -2.9664 - 2.8394i  -2.9664 - 2.8394i  -2.9664 - 2.8394i
```
```matlab
%% Product of two Toeplitz matrices yields a |TLMat| object
disp(TM1 * TM2)
```
```
  8x8 TLMat, displacement rank 4
```
```matlab
%% Convert |ToepMat| to full array
disp(full(TM1))
```
```
  Columns 1 through 5

   0.9947 + 0.3872i   1.7919 + 1.2852i  -0.0768 - 0.3662i  -0.2414 + 0.1604i  -2.0367 - 0.4784i
   0.1401 - 0.7831i   0.9947 + 0.3872i   1.7919 + 1.2852i  -0.0768 - 0.3662i  -0.2414 + 0.1604i
  -2.2595 - 0.8750i   0.1401 - 0.7831i   0.9947 + 0.3872i   1.7919 + 1.2852i  -0.0768 - 0.3662i
   1.6197 + 0.1019i  -2.2595 - 0.8750i   0.1401 - 0.7831i   0.9947 + 0.3872i   1.7919 + 1.2852i
  -0.5290 - 1.2252i   1.6197 + 0.1019i  -2.2595 - 0.8750i   0.1401 - 0.7831i   0.9947 + 0.3872i
  -1.9547 - 0.6110i  -0.5290 - 1.2252i   1.6197 + 0.1019i  -2.2595 - 0.8750i   0.1401 - 0.7831i
   1.0861 - 0.1590i  -1.9547 - 0.6110i  -0.5290 - 1.2252i   1.6197 + 0.1019i  -2.2595 - 0.8750i
  -2.0636 + 0.3248i   1.0861 - 0.1590i  -1.9547 - 0.6110i  -0.5290 - 1.2252i   1.6197 + 0.1019i

  Columns 6 through 8

   0.2668 + 0.0122i  -0.7960 - 0.9941i  -0.7643 - 0.8987i
  -2.0367 - 0.4784i   0.2668 + 0.0122i  -0.7960 - 0.9941i
  -0.2414 + 0.1604i  -2.0367 - 0.4784i   0.2668 + 0.0122i
  -0.0768 - 0.3662i  -0.2414 + 0.1604i  -2.0367 - 0.4784i
   1.7919 + 1.2852i  -0.0768 - 0.3662i  -0.2414 + 0.1604i
   0.9947 + 0.3872i   1.7919 + 1.2852i  -0.0768 - 0.3662i
   0.1401 - 0.7831i   0.9947 + 0.3872i   1.7919 + 1.2852i
  -2.2595 - 0.8750i   0.1401 - 0.7831i   0.9947 + 0.3872i
```
```matlab
%% Evaluate deg-6 Taylor polynomial of expm(TM1)
p = 1./factorial(6:-1:0);
E = polyvalm(p, TM1);  % No "full" arithmetic here!
disp(E); % Result is a TLMat
```
```
  8x8 TLMat, displacement rank 8
```
```matlab
EE = polyvalm(p, full(TM1)); % Compare with dense computation
disp(norm(E - EE, 'fro') / norm(EE, 'fro'));
```
```
   1.8832e-15
```
```matlab
%% Toeplitz-inverse is a |TLMat| of d-rank 2
disp(inv(TM1));
```
```
  8x8 TLMat, displacement rank 2
```
```matlab
%% Shifted inverse, useful for partial fraction expansion
sigma = exp(3i/4 * pi);
disp(inv(TM1 - sigma * toepeye(8)));
  8x8 TLMat, displacement rank 2
%% Solve a linear system
b = TM1 * ones(8,1);
x = TM1 \ b;
disp(norm(x - ones(8,1)));
```
```
   1.5422e-15
```
```matlab
%% It's also quite fast, using the GKO algorithm
[c,r,T] = random_toeplitz(4096,4096);
TMbig = ToepMat(c,r);
b = randn(4096,1);
tic; x = TMbig \ b; toc
```
```
Elapsed time is 0.947514 seconds.
```
```matlab
%% Compare with full arithmetic
tic; xx = T\b; toc
```
```
Elapsed time is 2.697644 seconds.
```
```matlab
%% Rational function evaluation, here: Cayley transform
e2 = zeros(16,1);
e2(2) = 1;
TM = ToepMat(-e2, e2); % Skew-symm
rT = TM.ratevalm([1,1], [1,-1]); % This is a |TLMat|
disp(norm(full(rT) - ( (full(TM) - eye(16)) \ (full(TM) + eye(16) ) )));
```
```
   1.5322e-15
```
```matlab
%% Can also be computed with \
rT = (TM - toepeye(16)) \ (TM + toepeye(16)); % Again a |TLMat|
disp(norm(full(rT) - ( (full(TM) - eye(16)) \ (full(TM) + eye(16) ) )));
```
```
   1.5304e-15
```
```matlab
%% Determinant is supported through GKO/LU factorization
TM = ToepMat([1,2,3,4]);
disp(det(TM))
```
```
  -20.0000
```
```matlab
%% Construct A TLMat object from a generator
G = orth(randn(8,4));
B = orth(randn(8,4)) * diag([10,1,0.01, 1e-12]);
TL1 = TLMat(G, B);  % drank is 4
TL2 = TL1.truncate_tol(1e-8); % remove ranks relatively smaller than 1e-8
fprintf('Trunctated from %d to %d, generator difference: %.2e\n', ...
    drank(TL1), drank(TL2), gennorm(TL1 - TL2));
```
```
Trunctated from 4 to 3, generator difference: 1.00e-12
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

Many details on the scope and implementation of TLComp are subject of an upcoming publication [3].


References
----


[1] Kailath, T., & Sayed, A. H. (Eds.). (1999). Fast reliable algorithms for matrices with structure. Society for Industrial and Applied Mathematics (SIAM), Philadelphia, PA. https://doi.org/10.1137/1.9781611971354

[2] A. Aricò and G. Rodriguez.  A fast solver for linear systems with
displacement structure.  Numer. Algorithms, 55(4):529-556, 2010.  DOI:
https://doi.org/10.1007/s11075-010-9421-x

[3] TLComp paper, in preparation
