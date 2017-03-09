%OUT1 drsolve numerical experiment #1, output routine
%   solution of random complex structured systems

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Feb 27, 2010

clear
load test1

uselatex = 0;	% set to 1 to produce a LaTeX table

if uselatex
	labs = {'Vandermonde','Toeplitz','Toeplitz+Hankel','Cauchy-like'};
	fprintf('\n\\begin{table}[htb]\n\\centering\n')
	fprintf('\\begin{tabular}{lcc}\n')
	fprintf('& errors & exec.~time \\\\\n\\hline\n')
	for i = 1:size(ERR,1)
		fprintf('%s', labs{i})
		[man1 esp1] = manesp(ERR(i,1));
		[man2 esp2] = manesp(ERR(i,2));
		fprintf(' & $%3.1f\\cdot 10^{%d}/%3.1f\\cdot 10^{%d}$', ...
			man1, esp1, man2, esp2);
		fprintf(' & $%.2f/%.2f$', TIM(i,:));
		fprintf(' \\\\\n')
	end
	fprintf('\\hline\n')
	fprintf('\\end{tabular}\n\\caption{insert caption}\n')
	fprintf('\\label{tab:tab1}\\end{table}\n\n')
else
	fprintf('\n\nError and execution time in the solution of complex structured linear systems\n')
	fprintf('\t[drsolve result] / [Matlab backslash result]  -  n = %g\n\n', n)
	fprintf('           Vandermonde        Toeplitz       Toeplitz+Hankel     Cauchy-like\n')
	fprintf('errors:  %.1e/%.1e   %.1e/%.1e   %.1e/%.1e   %.1e/%.1e\n', ...
		ERR(1,:), ERR(2,:), ERR(3,:), ERR(4,:))
	fprintf('time  :  %7.2f/%7.2f   %7.2f/%7.2f   %7.2f/%7.2f   %7.2f/%7.2f \n\n', ...
		TIM(1,:), TIM(2,:), TIM(3,:), TIM(4,:))
end
