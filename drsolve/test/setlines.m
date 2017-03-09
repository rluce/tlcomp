function setlines(hh,strt)
%SETLINES set line styles.
%
%   setlines(hh) sets the styles of the lines whose handles
%   are contained in the vector hh.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Feb 26, 2010

if nargin<2,  strt = 1;  end	% select first linestyle to use
styles = {'-','--',': '};
markers = {'o','s','*','x','+'};
sty = strt;
mar = 1;
for k = 1:length(hh)
	set( hh(k), 'linestyle', styles{sty})
	if ~strcmp(get(hh,'marker'),'none')
		set( hh(k), 'marker', markers{mar})
	end
	sty = sty + 1;
	if sty > length(styles),  sty = 1;  mar = mar+1;  end
	if mar > length(markers),  mar = 1;  end
end

