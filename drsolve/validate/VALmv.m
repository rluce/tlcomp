% VALmv

value = 1e-12;

% FTIMES
fprintf('     test FTIMES: ... ')
E = 0;
for n=2:30
    %
    F  = ftimes(n,'N');
    Fi = ftimes(n,'A');
    E  = max(E,norm( F*Fi-eye(n) ));
    clear F Fi
    %
    z  = exp(100i*rand(1));
    F  = ftimes(n,'N',z);
    Fi = ftimes(n,'A',z);
    E  = max(E,norm( F*Fi-eye(n) ));
    clear z F Fi
    %
    for j=1:30
        x  = complex(rand(n,j),randn(n,j));
        Fx = ftimes(x,'N');
        X  = ftimes(Fx,'A');
        E  = max(E,norm( x-X ));
        clear x Fx X
        %
        z  = exp(100i*rand(1));
        x  = sparse(complex(rand(n,j),randn(n,j)));
        Fx = ftimes(x,'N',z);
        X  = ftimes(Fx,'A',z);
        E  = max(E,norm( x-X ));
        clear z x Fx X
    end
end
if abs(E)<value, s='ok'; else s='KO'; end
disp([s ' (' num2str(E) ')']);
clear n j E s

% CTIMES
if exist('dct','file') && exist('idct','file')
	fprintf('     test CTIMES: ... ')
	E = 0;
	for n=2:30
	    %
	    C  = ctimes(n,'N');
	    Ct = ctimes(n,'T');
	    E  = max(E,norm( C*Ct-eye(n) ));
	    clear C Ct
	    %
	    for j=1:30
		x  = complex(rand(n,j),randn(n,j));
		Cx = ctimes(x,'N');
		X  = ctimes(Cx,'A');
		E  = max(E,norm( x-X ));
		clear x Cx X
	     end
	end
	if abs(E)<value, s='ok'; else s='KO'; end
	disp([s ' (' num2str(E) ')']);
	clear n j E s
else
	fprintf('     test CTIMES: SKIPPED\n')
end

% STIMES
if exist('dst','file')
	fprintf('     test STIMES: ... ')
	E = 0;
	for n=2:30
	    %
	    S  = stimes(n);
	    E  = max(E,norm( S*S-eye(n) ));
	    clear S
	    %
	    for j=1:30
		x  = sparse(complex(rand(n,j),randn(n,j)));
		Sx = stimes(x);
		X  = stimes(Sx);
		E  = max(E,norm( x-X ));
		clear x Sx X
	     end
	end
	if abs(E)<value, s='ok'; else s='KO'; end
	disp([s ' (' num2str(E) ')']);
	clear n j E s
else
	fprintf('     test STIMES: SKIPPED\n')
end

% CLTIMES+CL2FULL
fprintf('     test CLTIMES & CL2FULL: ... ')
E = 0;
for n=2:10
    for m=2:10
        for r=2:4
            t = complex(rand(m,1),rand(m,1));
            s = complex(rand(n,1),rand(n,1));
            G = complex(rand(m,r),rand(m,r));
            H = complex(rand(n,r),rand(n,r));
            for j=1:4
                x = complex(rand(n,j),rand(n,j));
                E = max(E,norm( cl2full(G,H,t,s)*x-cltimes(G,H,t,s,x) ));
            end
            clear t s G H x
        end
    end
end
if abs(E)<value, s='ok'; else s='KO'; end
disp([s ' (' num2str(E) ')']);
clear n m r j E s

clear value
