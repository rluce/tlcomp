function [col, row, T, b] = lps_example3(n)
% [c,r,T] = lps_example1(n)
%
% Example 3 from Lee, Hong-Kui Pang, and Hai-Wei Sun. Shift-invert Arnoldi
% approximation to the Toeplitz matrix exponential. SIAM J. Sci. Comput.,
% 32(2):774-792, 2010
%
% Time steps for expm(t*T) used in their work are t=0.5 and t=1.0

% Parameters used in their experiment
xi_min = -2;
xi_max = 2;
Delta_xi = (xi_max - xi_min) / (n+1);
K = 100;
nu = 0.25;
r = 0.05;
lambda = 0.1;
mu = -0.9;
sigma = 0.45;

kappa = exp( mu + sigma^2/2) - 1;

% TODO use Matlab's built in Gaussian pdf functionality
theta = @(eta) exp( - (eta-mu).^2 / (2*sigma^2) ) / (sqrt(2 * pi) * sigma);

% Initial condition for European call option
omega = @(xi) max( K*exp(xi) - K, 0 );

% Tridiagonal matrix corresponding to central differences
mid = nu^2/Delta_xi^2;
off = (2*r - 2*lambda*kappa - nu^2)/(4*Delta_xi);
D = diag((mid/2 - off) * ones(n-1,1), -1 ) + ...
    diag( (-mid - r - lambda) * ones(n,1) ) + ...
    diag((mid/2 + off) * ones(n-1,1), 1);

% Integral part of the equation
Ic = theta( ( 0:-1:(1-n) )' * Delta_xi );
Ir = theta( ( 0:(n-1) )  * Delta_xi );
I = Delta_xi * toeplitz(Ic, Ir);

% Assemble total matrix and rhs
T = D + lambda * I;
col = T(:,1);
row = T(1,:);
b = omega(linspace(xi_min, xi_max, n));

end