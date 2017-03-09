%TEST6	drsolve numerical experiment #6
%   solution of a Cauchy-like system with almost multiple knots

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 24, 2010

maxsiz = 2^12;		% max size for backslash and complete pivoting
tauvet = [1e-16 1e-14 1e-12 1e-10 1e-8]'; % parameters for clmk.m

n = 260;		% size
pivet = [1 2 4 5];	% pivoting techniques for tsolve

dr = 5;			% displacement rank
mult = dr;		% knots multiplicity
kb = 1;
rep = 1;

ntests = numel(tauvet);
seed = 0;
rand( 'state', seed)
randn( 'state', seed)

war = warning('off','drsolve:clsolve:multipleS');

% methods are: tsolve (var. pivoting), tsolve (collaps. knots), Matlab backslash
nmethods = numel(pivet) + 1 + 1;
ERR = zeros(ntests,nmethods);
TIM = zeros(ntests,nmethods);
CON = zeros(ntests,1);

fprintf('\ntau = ')

for i = 1:ntests
	tau = tauvet(i);
	fprintf('%g ', tau)
	for run = 1:rep
		[G H t s tc sc] = clmk(n,tau,dr,mult);
		sol = ones(n,kb);
		normsol = norm(sol,inf);
		b = cltimes(G,H,t,s,sol);

		% clsolve, various pivoting
		for j = 1:numel(pivet)
			piv = pivet(j);
			if (piv == 5) && (n > maxsiz)
				TIM(i,j) = nan;
				ERR(i,j) = nan;
			else
				tic
				x = clsolve( G, H, t, s, b, piv);
				TIM(i,j) = TIM(i,j) + toc;
				ERR(i,j) = ERR(i,j) + norm(x-sol,inf)/normsol;
			end
		end
		pos = numel(pivet);

		% clsolve, collapsed knots
		piv = 2;
		tic
		x = clsolve( G, H, tc, sc, b, piv);
		TIM(i,pos+1) = TIM(i,pos+1) + toc;
		ERR(i,pos+1) = ERR(i,pos+1) + norm(x-sol,inf)/normsol;
		pos = pos + 1;

		% Matlab backslash
		if n > maxsiz
			TIM(i,pos+1) = nan;
			ERR(i,pos+1) = nan;
			CON(i) = nan;
		else
			C = cl2full( G, H, t, s);
			CON(i) = cond(C);
			tic
			x = C \ b;
			TIM(i,pos+1) = TIM(i,pos+1) + toc;
			ERR(i,pos+1) = ERR(i,pos+1) + norm(x-sol,inf)/normsol;
		end
		pos = pos + 1;

	end
end
fprintf('\n')

ERR = ERR/rep;
TIM = TIM/rep;

warning(war)

clear C
save test6 ntests n tauvet pivet rep tau mult dr kb ERR TIM CON

out6

