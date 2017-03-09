% DRSolve: Displacement Rank Solver
% Solution of linear systems with displacement stucture.
% Version 1.0
%
% Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
% Email: {arico,rodriguez}@unica.it
%
% Solvers.
%   clsolve  - Solution of a Cauchy-like linear system.
%   tsolve   - Solution of a Toeplitz linear system.
%   tlsolve  - Solution of a Toeplitz-like linear system.
%   thsolve  - Solution of a Toeplitz+Hankel linear system.
%   thlsolve - Solution of a Toeplitz+Hankel-like linear system.
%   vsolve   - Solution of a Vandermonde linear system.
%   vlsolve  - Solution of a Vandermonde-like linear system.
%
% Conversion routines.
%   t2cl     - Converts a Toeplitz matrix to Cauchy-like.
%   t2tl     - Converts a Toeplitz matrix to Toeplitz-like.
%   tl2cl    - Converts a Toeplitz-like matrix to Cauchy-like.
%   th2cl    - Converts a square Toeplitz+Hankel matrix to Cauchy-like.
%   thl2cl   - Converts a square Toeplitz+Hankel-like matrix to Cauchy-like.
%   v2cl     - Converts a square Vandermonde matrix to Cauchy-like.
%   vl2cl    - Converts a square Vandermonde-like matrix to Cauchy-like.
%
% Conversion routines.
%   ftimes   - Fourier matrix product.
%   ctimes   - Cosine matrix product.
%   stimes   - Sine matrix product.
%   ttimes   - Toeplitz matrix product.
%   cltimes  - Cauchy-like matrix product.
%   cl2full  - Construct a Cauchy-like matrix from its generators.
%   nroots1  - n-th roots of unity.

