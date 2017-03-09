%OUT4 drsolve numerical experiment #4, output routine
%   growth of generators: Bini-Boito example (bbgcd)

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Feb 17, 2010

clear
load test4
seeall = 0;

fprintf('\n\nSolution of the Bini-Boito example\n\ttau=%g\n\n', tau)

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
cut = 12;
h1 = semilogy( (0:n-cut)', norh(1:end-cut,pv));
set(gca, 'fontsize', 14)
title('norm(H)')
legend(sleg{pv},2)
setlines(h1,3)
set(h1(1),'linestyle','-.')

fprintf('        [partial]   [Gu]     [total]  [backslash]\n')
fprintf('errors:  %.1e   %.1e   %.1e   %.1e\n\n', ERR([pv end]))

