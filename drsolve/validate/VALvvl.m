% VALvvl

value = 1e-7;
E = 0;
fprintf('     test V2CL, VL2CL, VSOLVE, VLSOLVE: ... ')
for n=2:30
    for j=1:4
        k = 0:n-1;
        y = exp(1i*(2*k+randn(1,n))*pi/n);
	eta = sqrt(y(1));
        b = complex(rand(n,j),rand(n,j));
        if rand(1)>0.5, b=real(b); end
        %
        V = vander(y);
        x1= V\b;
        x2= vsolve(y,eta,b);
        %
        E = max(E,norm( x1-x2 ));
        clear k y b V x1 x2 eta
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for dr = 1:min(n,5);
            t = nroots1(n,1);
            for k=1:3
                G = complex(rand(n,dr),rand(n,dr));
                H = complex(rand(n,dr),rand(n,dr));
                sol=complex(rand(n,k),rand(n,k));
                eta=exp(100i*randn(1));
                s = conj(nroots1(n,eta));
                %
                switch mod(j,3)
                    case 0
                        G=real(G);H=real(H);sol=real(sol);
                    case 1
                    if rand(1)>0.5, G=real(G); end
                    if rand(1)>0.5, H=real(H); end
                    if rand(1)>0.5, sol=real(sol); end
                end
                b = cltimes(G,H,t,s,sol);
                % C2VL
                H = ftimes(H,'N',eta);
                x = vlsolve(t,eta,G,H,b,1);
                % VL2C
                x = ftimes(x,'A',eta);
                E = max(E,norm( x-sol ));
                clear G H sol eta s b x
            end
        end
    end
end
clear dr k t

if abs(E)<value, s='ok'; else s='KO'; end
disp([s ' (' num2str(E) ')']);

clear E n j s value
