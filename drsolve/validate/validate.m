% validate

rand('seed',0)
randn('seed',0)
disp(' ')
disp('Validation script for the DRSOLVE package.')
disp(' ')

switch exist('clsolve','file')
    case 0
        disp('- DRSOLVE is not properly installed. CLSOLVE is not in your path!');
        return
    case 2
        disp('- DRSOLVE is installed, but MEX file are missing. Try to run installmex.')
        disp('- Computation may be slow.')
        fprintf('  press a key to continue...'), pause, fprintf('\n')
    case 3
        disp('- DRSOLVE is properly installed.')
    otherwise
        disp('Unable to check. STOP')
        return
end

if exist('dct','file') && exist('idct','file') && exist('dst','file')
    disp('- DCT, IDCT and DST are available on your system.')
else
    disp('- DCT, IDCT and/or DST are not available on your system.')
    disp('  They are needed for THSOLVE and THLSOLVE only.')
    fprintf('  press a key to continue...'), pause, fprintf('\n')
end

disp('- now some tests will be run to check that the routines work...')

disp('  1) fast products & co.')
VALmv

war = warning('off','drsolve:clsolve:missingMEX');

disp('  2) Toeplitz and Toeplitz-like computation')
VALttl

if exist('dct','file') && exist('idct','file') && exist('dst','file')
    disp('  3) Toeplitz+Hankel and Toeplitz+Hankel-like computation')
    VALththl
else
    disp('  3) Toeplitz+Hankel and Toeplitz+Hankel-like computation SKIPPED')
end

disp('  4) Vandermonde and Vandermonde-like computation')
VALvvl

warning(war)
clear war

