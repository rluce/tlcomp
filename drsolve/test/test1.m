%TEST1	drsolve numerical experiment #1
%   solution of random complex structured systems

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

l2nvet = [11];
nvet = 2.^l2nvet;
pivet = [1];

rep = 1;
dr = 5;

ntests = numel(nvet);
seed = 0;
rand( 'state', seed)
randn( 'state', seed)

ERR = zeros(4,2);
TIM = zeros(4,2);
%CND = zeros(4,1);

fprintf('\nn     = %6d\n', nvet(end))
fprintf('rep   = %6d\n',rep)

n = nvet(1);
piv = pivet(1);
sol = ones(n,1);
normsol = norm(sol,inf);

for run = 1:rep
	fprintf('test %g: ',run)

	% Vandermonde
	fprintf('V, ');
	w = exp((2*(0:n-1)'+1+randn(n,1))*1i*pi/n);
	eta = 1;
	A = vander(w);
	%CND(1) = CND(1)+cond(A);
	b = A*sol;
	tic
	x = vsolve( w, eta, b, piv);
	TIM(1,1) = TIM(1,1) + toc;
	ERR(1,1) = ERR(1,1) + norm(x-sol,inf)/normsol;
	tic
	x = A \ b;
	TIM(1,2) = TIM(1,2) + toc;
	ERR(1,2) = ERR(1,2) + norm(x-sol,inf)/normsol;

	% Toeplitz
	fprintf('T, ');
	c = rand(n,1) + 1i*rand(n,1);
	r = rand(n,1) + 1i*rand(n,1);
	r(1) = c(1);
	A = toeplitz(c,r);
	%CND(2) = CND(2)+cond(A);
	b = A*sol;
	tic
	x = tsolve( c, r, b, piv);
	TIM(2,1) = TIM(2,1) + toc;
	ERR(2,1) = ERR(2,1) + norm(x-sol,inf)/normsol;
	tic
	x = A \ b;
	TIM(2,2) = TIM(2,2) + toc;
	ERR(2,2) = ERR(2,2) + norm(x-sol,inf)/normsol;

	% Toeplitz+Hankel
	fprintf('TH, ');
	h1 = rand(n,1) + 1i*rand(n,1);
	h2 = rand(n,1) + 1i*rand(n,1);
	h2(1) = h1(n);
	A = toeplitz(c,r)+hankel(h1,h2);
	%CND(3) = CND(3)+cond(A);
	b = A*sol;
	tic
	x = thsolve( c, r, h1, h2, b, piv);
	TIM(3,1) = TIM(3,1) + toc;
	ERR(3,1) = ERR(3,1) + norm(x-sol,inf)/normsol;
	tic
	x = A \ b;
	TIM(3,2) = TIM(3,2) + toc;
	ERR(3,2) = ERR(3,2) + norm(x-sol,inf)/normsol;

	% Cauchy
	fprintf('CL\n');
	G = randn(n,dr) + 1i*randn(n,dr);
	H = randn(n,dr) + 1i*randn(n,dr);
	t = nroots1(n,1);
	s = nroots1(n,-1);
	A = cl2full(G,H,t,s);
	%CND(4) = CND(4)+cond(A);
	b = A*sol;
	tic
	x = clsolve( G, H, t, s, b, piv);
	TIM(4,1) = TIM(4,1) + toc;
	ERR(4,1) = ERR(4,1) + norm(x-sol,inf)/normsol;
	tic
	x = A \ b;
	TIM(4,2) = TIM(4,2) + toc;
	ERR(4,2) = ERR(4,2) + norm(x-sol,inf)/normsol;

end
fprintf('\n')

ERR = ERR/rep;
TIM = TIM/rep;
%CND = CND/rep;

clear A
save test1 n rep dr piv ERR TIM %CND

out1


