%OUT3 drsolve numerical experiment #3, output routine
%   solution of Gaussian real Toeplitz systems

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Feb 11, 2010

clear
load test3
seeall = 0;

fprintf('\n\nError and execution time in the solution of Gaussian real Toeplitz linear systems\n\tsigma=%g\n\n', sigma)

sleg{1}='no piv.';
sleg{2}='partial';
sleg{3}='mult. s';
sleg{4}='S&B';
sleg{5}='Gu';
sleg{6}='total';

if seeall
figure(10)
clf
h10 = loglog( nvet, ERR(:,:), '-o');
set(gca, 'fontsize', 14)
title('All together!')
xlim( nvet([1 ntests]))
set( gca, 'xtick', nvet)
legend(sleg{pivet+1},'toms729(2)','toms729(10)','T+H partial','T no MEX','LU',4)
setlines(h10)
end

figure(1)
clf
h1 = loglog( nvet, ERR(:,[2 5 9 8 11]), '-o');
set(gca, 'fontsize', 14)
title('Relative error')
xlim( nvet([1 ntests]))
set( gca, 'xtick', nvet)
legend('T (partial)','T (Gu)','T+H (partial)','toms729','\it{backslash}',2)
setlines(h1)

figure(2)
clf
alfan2 = nvet'.^2;
alfan2 = alfan2 / alfan2(end) * max(max(TIM)) / 10^(2);
h2 = loglog( nvet, [TIM(:,[2 9 8 10 11]) alfan2], '-o');
set(gca, 'fontsize', 14)
title('Execution time')
xlim( nvet([1 ntests]))
set( gca, 'xtick', nvet)
legend('T (partial)','T+H (partial)','toms729','T (no MEX)','\it{backslash}','\alpha n^2',4)
setlines(h2)

figure(3)
clf
h3 = loglog( nvet, ERR(:,[1 2 4 5 6]), '-o');
set(gca, 'fontsize', 14)
title('Relative error')
xlim( nvet([1 ntests]))
set( gca, 'xtick', nvet)
legend(sleg{[1 2 4 5 6]},2)
setlines(h3)

