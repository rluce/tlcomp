%OUT6 drsolve numerical experiment #6, output routine
%   solution of a Cauchy-like system with almost multiple knots

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Feb 18, 2010

clear
load test6
seeall = 0;

fprintf('\n\nError in the solution of a Cauchy-like system with almost multiple knots\n\n')

sleg{1}='no piv.';
sleg{2}='partial';
sleg{3}='mult. s';
sleg{4}='S&B';
sleg{5}='Gu';
sleg{6}='total';

if seeall
figure(10)
clf
h10 = loglog( tauvet, ERR(:,:), '-o');
set(gca, 'fontsize', 14)
title('All together!')
legend(sleg{pivet+1},'coll. knots','\it{backslash}',1)
setlines(h10)
end

figure(1)
clf
pv = [2 5 6];
h1 = loglog( tauvet, ERR(:,[1 3:6]), '-o');
set(gca, 'fontsize', 14)
title('Relative error')
xlabel('\tau values')
legend(sleg{[2 5 6]},'coll. knots','\it{backslash}',1)
setlines(h1)
