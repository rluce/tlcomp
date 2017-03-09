function installmex(forceversion)
%INSTALLMEX Install the appropriate MEX version of clsolve, if available

%   installmex('x.y') will force installation for version x.y of Matlab

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Oct 20, 2010

clear functions
comp = computer;
verstr = ver('matlab');
filnam = '';

if nargin > 0
	verstr.Version = forceversion;
end

switch comp
    case 'PCWIN'
        switch verstr.Version
            case '7.4'
                filnam = '../src/win32-matlab74/clsolve.mexw32';
            case '7.7'
                filnam = '../src/win32-matlab77/clsolve.mexw32';
            case {'7.9','7.10'}
                filnam = '../src/win32-matlab79/clsolve.mexw32';
        end
    case 'GLNX86'
        filnam = '../src/linux32-matlab/clsolve.mexglx';
    case 'GLNXA64'
        switch verstr.Version
            case {'7.4','7.5','7.6','7.7'}
                filnam = '../src/linux64-matlab74/clsolve.mexa64';
            case {'7.8','7.9'}
                filnam = '../src/linux64-matlab78/clsolve.mexa64';
            case {'7.10','7.11'}
                filnam = '../src/linux64-matlab710/clsolve.mexa64';
            
        end
    case 'MACI'
        switch verstr.Version
            case '7.7'
                filnam = '../src/macx86-matlab77/clsolve.mexmaci';
        end
end

if isempty(filnam)
    disp('The MEX version of clsolve is not available for this OS/Matlab version.')
else
    if ~copyfile( filnam, '../');
        error('File not found.')
    end
end

