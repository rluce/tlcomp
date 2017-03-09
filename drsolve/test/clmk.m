function [G H t s tc sc] = clmk(n,tau,dr,mult,graf)
%CLMK	Cauchy-like matrix with "almost" multiple knots
%   [G,H,t,s,tc,sc]=clmk(n,tau,dr,mult) returns the displacement structure 
%   of a square Cauchy-like matrix of displacement rank dr and dimension n 
%   whose knots are repeated mult time with a random perturbation scaled 
%   by tau.  The vectors tc and sc optionally return the knots without
%   perturbation.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it

%   Last revised Jun 17, 2010

if nargin<2 || isempty(tau),  tau = 1e-8;  end
if nargin<3 || isempty(dr),  dr = 4;  end
if nargin<4 || isempty(mult),  mult = dr;  end
if nargin<5 || isempty(graf),  graf = 0;  end

if mult>dr,  error('singular matrix'),  end

%rand('seed',0)
%randn('seed',0)
G = rand(n,dr);
H = rand(n,dr);
m = ceil(n/mult);
t = zeros(m*mult,1);
s = zeros(m*mult,1);
t(1:m) = nroots1(m,1);
s(1:m) = nroots1(m,-1);
for i = 2:mult
	t((i-1)*m+1:i*m) = t(1:m).*(1+tau*complex(randn(m,1),randn(m,1)));
	s((i-1)*m+1:i*m) = s(1:m).*(1+tau*complex(randn(m,1),randn(m,1)));
end
t = t(1:n);
s = s(1:n);
tc = repmat(t(1:m),mult,1);
sc = repmat(s(1:m),mult,1);
tc = tc(1:n);
sc = sc(1:n);
if graf
	plot( s, 'o')
	set(gca, 'fontsize', 14)
	axis square
	axis([-1.2 1.2 -1.2 1.2])
end

