% VALththl

value = 1e-7;

fprintf('     test TH2CL, THL2CL, THSOLVE, THLSOLVE: ... ')
warn = warning;
warning off
E = 0;
for n=4:20
    for j=1:6
        % test THSOLVE
        for k=1:4
            rT = complex(rand(n,1),rand(n,1));
            cT = complex(rand(n,1),rand(n,1));
            rH = complex(rand(n,1),rand(n,1));
            cH = complex(rand(n,1),rand(n,1));
            b  = complex(rand(n,k),rand(n,k));
            switch mod(j,3)
                case 0
                    rT=real(rT);cT=real(cT);rH=real(rH);cH=real(cH);b=real(b);
                case 1
                    if rand(1)>0.5, rT=real(rT); end
                    if rand(1)>0.5, cT=real(cT); end
                    if rand(1)>0.5, rH=real(rH); end
                    if rand(1)>0.5, cH=real(cH); end
                    if rand(1)>0.5, b =real(b);  end
            end            
            T = toeplitz(cT,rT);
            H = hankel(cH,rH);
            x1= (T+H)\b;
            x2= thsolve(cT,rT,cH,rH,b);
            E = max(E,norm( x1-x2 ));
            clear rT cT rH cH b T H x1 x2
        end
        % test THLSOLVE
        t = 2*cos((1:n).'*(pi/(n+1)));
        s = 2*cos((0:n-1).'*(pi/n));
        for dr = 1:min(n,5);
            for k=1:3
                G = complex(rand(n,dr),rand(n,dr));
                H = complex(rand(n,dr),rand(n,dr));
                sol=complex(rand(n,k),rand(n,k));
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
                % C2(T+H)L
                G = stimes(G);
                H = ctimes(H,'N');
                b = stimes(b);
                % (T+H)L2C
                x = thlsolve(G,H,b,1);
                x = ctimes(x,'T');
                E = max(E,norm( x-sol ));
                clear G H sol b x
            end
        end
        clear t s
    end
end

if abs(E)<value, s='ok'; else s='KO'; end
disp([s ' (' num2str(E) ')']);
warning(warn)

clear n j k dr E warn value s
