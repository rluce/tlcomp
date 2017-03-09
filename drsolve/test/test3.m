%TEST3	drsolve numerical experiment #3
%   solution of Gaussian real Toeplitz systems

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 19, 2010

l2nvet = [7:12];	% range of exponents (size is 2^k)
maxsiz = 2^10;		% max size for backslash and complete pivoting
%l2nvet = [8:15];	% values used in the paper
%maxsiz = 2^12;		% values used in the paper

nvet = 2.^l2nvet;	% size is 2^k
pivet = [0 1 2 3 4 5];	% pivoting techniques for tsolve
pmaxv = [2 10];		% pmax for toms729
sigma = .3;		% parameter for the Gaussian matrix

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
TIM = zeros(ntests,nmethods);

fprintf('\nmax n = %d\n\n', nvet(end))

for i = 1:ntests
	n = nvet(i);
	fprintf('%d ', n)
	for run = 1:rep
		c = sqrt(sigma/2/pi) * exp(-sigma/2*(0:n-1)'.^2);
		r = c;
		r(1) = c(1);
		sol = ones(n,kb);
		normsol = norm(sol,inf);
		b = ttimes(c,r,sol);

		% tsolve, various pivoting
		for j = 1:numel(pivet)
			piv = pivet(j);
			if (piv == 5) && (n > maxsiz)
				TIM(i,j) = nan;
				ERR(i,j) = nan;
			else
				tic
				x = tsolve( c, r, b, piv);
				TIM(i,j) = TIM(i,j) + toc;
				ERR(i,j) = ERR(i,j) + norm(x-sol,inf)/normsol;
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
			else
				TIM(i,pos+p) = nan;
				ERR(i,pos+p) = nan;
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
		pos = pos + 1;

		% tsolve, partial pivoting, no MEX
		piv = 1;
		if n > maxsiz
			TIM(i,pos+1) = nan;
			ERR(i,pos+1) = nan;
		else
			tic
			x = tsolvenomex( c, r, b, piv);
			TIM(i,pos+1) = TIM(i,pos+1) + toc;
			ERR(i,pos+1) = ERR(i,pos+1) + norm(x-sol,inf)/normsol;
		end
		pos = pos + 1;

		% Matlab backslash
		if n > maxsiz
			TIM(i,pos+1) = nan;
			ERR(i,pos+1) = nan;
		else
			T = toeplitz(c,r);
			tic
			x = T \ b;
			TIM(i,pos+1) = TIM(i,pos+1) + toc;
			ERR(i,pos+1) = ERR(i,pos+1) + norm(x-sol,inf)/normsol;
		end
		pos = pos + 1;

	end
end
fprintf('\n')

ERR = ERR/rep;
TIM = TIM/rep;

clear T
save test3 ntests sigma l2nvet nvet pivet pmaxv rep dr kb ERR TIM

out3

