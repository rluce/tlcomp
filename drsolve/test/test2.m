%TEST2	drsolve numerical experiment #2
%   solution of random real Toeplitz systems

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Jul 27, 2010

l2nvet = [7:12];	% range of exponents (size is 2^k)
maxsiz = 2^10;		% max size for backslash and complete pivoting
%l2nvet = [8:15];	% values used in the paper
%maxsiz = 2^12;		% values used in the paper

nvet = 2.^l2nvet;	% size is 2^k
pivet = [0 1 2 3 4 5];	% pivoting techniques for tsolve
pmaxv = [2 10];		% pmax for toms729

dr = 2;
kb = 1;
rep = 1;

ntests = numel(nvet);
seed = 0;
rand( 'state', seed)
randn( 'state', seed)

% check if toms729 is available
havetoms = (exist('toms729','file') == 3);
if ~havetoms, disp('TOMS729 is missing, test will be skipped'),  end

% methods are: tsolve (var. pivoting), toms729 (var. pmax), 
%              thsolve (partial pivoting), tsolve (no mex), Matlab backslash
nmethods = numel(pivet) + numel(pmaxv) + 1 + 1 + 1;
ERR = zeros(ntests,nmethods);
RES = zeros(ntests,nmethods);
TIM = zeros(ntests,nmethods);

fprintf('\nmax n = %d\n\n', nvet(end))

for i = 1:ntests
	n = nvet(i);
	fprintf('%d ', n)
	for run = 1:rep
		c = rand(n,1);
		r = rand(n,1);
		r(1) = c(1);
		sol = ones(n,kb);
		normsol = norm(sol,inf);
		b = ttimes(c,r,sol);
		normb = norm(b,inf);

		% tsolve, various pivoting
		for j = 1:numel(pivet)
			piv = pivet(j);
			if (piv == 5) && (n > maxsiz)
				TIM(i,j) = nan;
				ERR(i,j) = nan;
				RES(i,j) = nan;
			else
				tic
				x = tsolve( c, r, b, piv);
				TIM(i,j) = TIM(i,j) + toc;
				ERR(i,j) = ERR(i,j) + norm(x-sol,inf)/normsol;
				RES(i,j) = RES(i,j) + norm(b-ttimes(c,r,x),inf)/normb;
			end
		end
		pos = numel(pivet);

		% toms729, various pmax
		for p = 1:numel(pmaxv)
			pmax = pmaxv(p);
			if havetoms
				tic
				x = toms729( c, r, b, pmax);
				TIM(i,pos+p) = TIM(i,pos+p) + toc;
				ERR(i,pos+p) = ERR(i,pos+p) + norm(x-sol,inf)/normsol;
				RES(i,pos+p) = RES(i,pos+p) + norm(b-ttimes(c,r,x),inf)/normb;
			else
				TIM(i,pos+p) = nan;
				ERR(i,pos+p) = nan;
				RES(i,pos+p) = nan;
			end
		end
		pos = pos + numel(pmaxv);

		% thsolve, partial pivoting, zero Hankel part
		cH = zeros(n,1);
		rH = zeros(n,1);
		piv = 1;
		tic
		x = thsolve( c, r, cH, rH, b, piv);
		TIM(i,pos+1) = TIM(i,pos+1) + toc;
		ERR(i,pos+1) = ERR(i,pos+1) + norm(x-sol,inf)/normsol;
		RES(i,pos+1) = RES(i,pos+1) + norm(b-ttimes(c,r,x),inf)/normb;
		pos = pos + 1;

		% tsolve, partial pivoting, no MEX
		piv = 1;
		if n > maxsiz
			TIM(i,pos+1) = nan;
			ERR(i,pos+1) = nan;
			RES(i,pos+1) = nan;
		else
			tic
			x = tsolvenomex( c, r, b, piv);
			TIM(i,pos+1) = TIM(i,pos+1) + toc;
			ERR(i,pos+1) = ERR(i,pos+1) + norm(x-sol,inf)/normsol;
			RES(i,pos+1) = RES(i,pos+1) + norm(b-ttimes(c,r,x),inf)/normb;
		end
		pos = pos + 1;

		% Matlab backslash
		if n > maxsiz
			TIM(i,pos+1) = nan;
			ERR(i,pos+1) = nan;
			RES(i,pos+1) = nan;
		else
			T = toeplitz(c,r);
			tic
			x = T \ b;
			TIM(i,pos+1) = TIM(i,pos+1) + toc;
			ERR(i,pos+1) = ERR(i,pos+1) + norm(x-sol,inf)/normsol;
			RES(i,pos+1) = RES(i,pos+1) + norm(b-ttimes(c,r,x),inf)/normb;
		end
		pos = pos + 1;

	end
end
fprintf('\n')

ERR = ERR/rep;
RES = RES/rep;
TIM = TIM/rep;

clear T
save test2 ntests l2nvet nvet pivet pmaxv rep dr kb ERR RES TIM

out2

