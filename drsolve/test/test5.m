%TEST5	drsolve numerical experiment #5
%   growth of generators: Sweet-Brent example (sbbad)

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Feb 18, 2010

n = 128;		% size of problem

pivet = [0:5];
tau = 1e-12;		% parameter for sbbad
gurep = 10;		% frequency of Gu's QR step

dr = 2;
kb = 1;

seed = 0;
rand( 'state', seed)
randn( 'state', seed)

ERR = zeros(2,numel(pivet)+1);
norg = zeros(n+1,numel(pivet));
norh = zeros(n+1,numel(pivet));

fprintf('\nn = %d\n\n', n)

[G H t s a f] = sbbad(n,tau);
sol = ones(n,1);
C = cl2full(-f,a,t,s);		% use true displacement for accuracy
b = C*sol;

for j = 1:numel(pivet)
	piv = pivet(j)+gurep*10;
	[x rc p q norg(:,j) norh(:,j)] = clsolvenorm( G, H, t, s, b, piv);
	ERR(1,j) = norm(x-sol,inf)/norm(sol,inf);
	x = clsolvenorm( -f, a, t, s, b, piv);
	ERR(2,j) = norm(x-sol,inf)/norm(sol,inf);
end

C1 = cl2full(G,H,t,s);
[L U P] = lu(C1);
x = U \ (L \ (P*b));
ERR(1,end) = norm(x-sol,inf)/norm(sol,inf);
[L U P] = lu(C);
x = U \ (L \ (P*b));
ERR(2,end) = norm(x-sol,inf)/norm(sol,inf);

clear C C1 L U P
save test5 n pivet tau gurep dr kb norg norh ERR

out5

