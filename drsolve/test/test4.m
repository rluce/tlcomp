%TEST4	drsolve numerical experiment #4
%   growth of generators: Bini-Boito example (bbgcd)

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Feb 12, 2010

n = 512;		% size of problem

pivet = [0:5];
tau = 1e-10;		% parameter for bbgcd
gurep = 10;		% frequency of Gu's QR step

dr = 2;
kb = 1;

seed = 0;
rand( 'state', seed)
randn( 'state', seed)

ERR = zeros(1,numel(pivet)+1);
norg = zeros(n+1,numel(pivet));
norh = zeros(n+1,numel(pivet));

fprintf('\nn = %d\n\n', n)

[G H t s] = bbgcd(n,20,tau);
sol = ones(n,1);
C = cl2full(G,H,t,s);
b = C*sol;

for j = 1:numel(pivet)
	piv = pivet(j)+gurep*10;
	[x rc p q norg(:,j) norh(:,j)] = clsolvenorm( G, H, t, s, b, piv);
	ERR(j) = norm(x-sol,inf)/norm(sol,inf);
end

[L U P] = lu(C);
x = U \ (L \ (P*b));
ERR(end) = norm(x-sol,inf)/norm(sol,inf);

clear C L U P
save test4 n pivet tau gurep dr kb norg norh ERR

out4

