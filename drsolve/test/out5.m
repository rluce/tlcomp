%OUT5 drsolve numerical experiment #5, output routine
%   growth of generators: Sweet-Brent example (sbbad)

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 23, 2010

clear
load test5
seeall = 0;

uselatex = 0;	% set to 1 to produce a LaTeX table

fprintf('\n\nSolution of the Sweet-Brent example\n\ttau=%g\n\n', tau)

sleg{1}='no piv.';
sleg{2}='partial';
sleg{3}='mult. s';
sleg{4}='S&B';
sleg{5}='Gu';
sleg{6}='total';

if seeall
figure(10)
clf
subplot(1,2,1)
h101 = semilogy( (0:n)', norg);
set(gca, 'fontsize', 14)
xlim( [0 n])
yl1 = ylim;
title('norm(G)')
setlines(h101)
subplot(1,2,2)
h102 = semilogy( (0:n)', norh);
set(gca, 'fontsize', 14)
title('norm(H)')
xlim( [0 n])
yl2 = ylim;
legend(sleg{pivet+1},3)
setlines(h102)
yl = [min([yl1(1) yl2(1)]) max([yl1(2) yl2(2)])];
subplot(1,2,1)
ylim( yl)
subplot(1,2,2)
ylim( yl)
end

figure(1)
clf
pv = [2 5 6];
cut = 8;
subplot(1,2,1)
h11 = semilogy( (0:n-cut)', norg(1:end-cut,pv));
set(gca, 'fontsize', 14)
xlim([0 n-cut])
title('norm(G)')
legend(sleg{pv},4)
setlines(h11,3)
set(h11(1),'linestyle','-.')
subplot(1,2,2)
h21 = semilogy( (0:n-cut)', norh(1:end-cut,pv));
set(gca, 'fontsize', 14)
xlim([0 n-cut])
title('norm(H)')
setlines(h21,3)
set(h21(1),'linestyle','-.')

if uselatex
	labs = {'partial','Gu','total','\emph{backslash}'};
	fprintf('\n\\begin{table}[htb]\n\\centering\n')
	fprintf('\\begin{tabular}{lcccc}\n')
	fprintf('& partial & Gu & total & \\emph{backslash} \\\\\n\\hline\n')
	fprintf('S-B displacement')
	for i = 1:size(pv,2)
		[man1 esp1] = manesp(ERR(1,pv(i)));
		fprintf(' & $%3.1f\\cdot 10^{%d}$', man1, esp1);
	end
	[man1 esp1] = manesp(ERR(1,end));
	fprintf(' & $%3.1f\\cdot 10^{%d}$ \\\\\n', man1, esp1);
	fprintf('alternative displ.')
	for i = 1:size(pv,2)
		[man1 esp1] = manesp(ERR(2,pv(i)));
		fprintf(' & $%3.1f\\cdot 10^{%d}$', man1, esp1);
	end
	[man1 esp1] = manesp(ERR(2,end));
	fprintf(' & $%3.1f\\cdot 10^{%d}$ \\\\\n', man1, esp1);
	fprintf('\\hline\n')
	fprintf('$\\max|G^{(k)}|/|G^{(0)}|$')
	for i = 1:size(pv,2)
		fprintf(' & $%4.1f$', max(norg(:,pv(i))/norg(1,pv(i))));
	end
	fprintf(' \\\\\n');
	fprintf('$\\max|H^{(k)}|/|H^{(0)}|$')
	for i = 1:size(pv,2)
		fprintf(' & $%4.1f$', max(norh(:,pv(i))/norh(1,pv(i))));
	end
	fprintf(' \\\\\n');
	fprintf('\\hline\n')
	fprintf('\\end{tabular}\n\\caption{insert caption}\n')
	fprintf('\\label{tab:tab5}\\end{table}\n\n')
else
	fprintf('errors        [partial]   [Gu]     [total]  [backslash]\n')
	fprintf('S-B displ.  :  %.1e   %.1e   %.1e   %.1e\n', ERR(1,[pv end]))
	fprintf('alt. displ. :  %.1e   %.1e   %.1e   %.1e\n\n', ERR(2,[pv end]))
	fprintf('max(norm(G)):  %5.1f      %5.1f    %5.1f\n', max(norg(:,pv))./norg(1,pv))
	fprintf('max(norm(H)):  %5.1f      %5.1f    %5.1f\n\n', max(norh(:,pv))./norh(1,pv))
end

