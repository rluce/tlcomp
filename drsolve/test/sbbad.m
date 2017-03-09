function [G,H,t,s,a,f] = sbbad(n,tau)
%SBBAD	Sweet&Brent Cauchy-like example
%   [G,H,t,s]=sbbad(n,r) returns the displacement structure of a square
%   Cauchy-like matrix of displacement rank 2 and dimension n, described 
%   in the paper
%   D.R. Sweet, R.P. Brent, Error analysis of a fast partial pivoting 
%   method for structured matrices. In Advanced Signal Processing Algorithms, 
%   F.T. Luk, Ed. Vol. 2563, SPIE, San Diego, 266--280 (1995).

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it

%   Last revised Feb 17, 2010

%a = randn(n,1);
%f = randn(n,1);
a = ones(n,1);
f = (-1).^(1:n)';
a = a/norm(a);
f = f/norm(f);
f = f*tau;
G = [a a+f];
H = [a -a];
t = nroots1(n,1);
s = nroots1(n,-1);

