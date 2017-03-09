function [G,H,t,s] = bbgcd(ndim,rankdef,tau)
%BBGCD	"almost" rank deficient Cauchy-like matrix
%   [G,H,t,s]=bbgcd(n,r) returns the displacement structure of a square
%   Cauchy-like matrix of displacement rank 2, dimension n and rank n-r
%   (r<=n/2).
%   [G,H,t,s]=bbgcd(n,r,tau) returns a full rank matrix by adding to the
%   generators G and H a perturbation scaled by tau.
%
%   This matrix exhibits a huge generators growth during the application of the
%   generalized Schur algorithm.
%
%   The example is taken from the paper: D.A. Bini, P. Boito, A fast algorithm
%   for approximate polynomial gcd based on structured matrix computations
%   (2008), submitted.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%   Thanks to Paola Boito for providing the code.
%
%   Last revised Feb 17, 2010

if rankdef > ndim/2,  error('This routine requires r<=n/2.'),  end

%--
d3=rankdef+1;
d1=floor(ndim/2)-d3+2;
d2=ceil(ndim/2)-d3+2;
%--
%randn( 'state', 0)
p=randn(1,d1);
q=randn(1,d2);
r=randn(1,d3);
%--
f=conv(p,r);
g=conv(q,r);
%--
n = d1+d3-2;
m = d2+d3-2;
N = n+m;
%--
G2=zeros(2,N);
G2(1,N-m:N)=g;
G2(1,1:n)=G2(1,1:n)-f(2:n+1);
G2(1,N)=G2(1,N)+f(1);
G2(2,N-n:N)=f;
G2(2,1:m)=G2(2,1:m)-g(2:m+1);
G2(2,N)=G2(2,N)+g(1);
G1=zeros(N,2);
G1([1;m+1],:)=eye(2);
%--
d0=exp(1i*pi/N*(0:N-1));
D1=d0.^2;
D2=exp(1i*pi/N)*D1;
sN=sqrt(N);
FD0=sN*ifft(diag(d0));
G1_c=sN*ifft(G1);
G2_c=G2*FD0';
G = G1_c+tau*randn(ndim,2);
H = G2_c'+tau*randn(ndim,2);
t = D1.';
s = D2.';

