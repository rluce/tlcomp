% VALttl

value = 1e-8;

fprintf('     test T2TL, TL2CL, T2CL, TSOLVE, TLSOLVE: ... ')
warn = warning;
warning off
E = 0;
for m=2:20
    for n=2:20
        r = complex(rand(1,n),rand(1,n));
        c = complex(rand(m,1),rand(m,1));
        %
        if rand(1)>0.5, r = real(r); end
        if rand(1)>0.5, c = real(c); end
        if rand(1)>0.5, r(1) = c(1); end
        %
        T = toeplitz(c,r);
        [G,H,t,s,xi,eta] = t2cl(c,r);
        C = cl2full(G,H,t,s);
        % test CL2FULL
        E = max(E,norm( diag(t)*C-C*diag(s) - G*H' ));
        % test T2CL !
        E = max(E,norm( C - ftimes( ftimes(T,'A',xi)','A',eta)' ));
        clear C
        if m==n
            for d=1:4
                % test TSOLVE
                B = complex(rand(n,d),rand(n,d));
                if rand(1)>0.5, B=real(B); end
                x1= tsolve(c,r,B);
                x2= T\B;
                E = max(E,norm( x1-x2 ));
            end
        end
        % T2TL+TL2CL+T2CL
        % 1
        [G1,H1,xi1,eta1] = t2tl(c,r);
        E = max([E,abs(xi-xi1),abs(eta-eta1)]);
        Z_xi  = sparse([2:m,1],1:m,[ones(1,m-1),xi ],m,m);
        Z_eta = sparse([2:n,1],1:n,[ones(1,n-1),eta],n,n);
        E = max(E,norm(Z_xi*T-T*Z_eta-G1*H1'));
        % 2
        [G2,H2,t2,s2] = tl2cl(G1,H1,xi1,eta1);
        E = max([E,norm(G-G2),norm(H-H2),norm(t-t2),norm(s-s2)]);
        % ancora test..
        xi  = exp(100i*rand(1));
        eta = exp(100i*rand(1));
        [G3,H3,t3,s3] = t2cl(c,r,xi,eta);
        [G4,H4] = t2tl(c,r,xi,eta); %!
        [G5,H5,t5,s5] = tl2cl(G4,H4,xi,eta);
        E = max([E,norm(G3-G5),norm(H3-H5),norm(t3-t5),norm(s3-s5)]);
        if m==n
            for d=1:4
                % test TLSOLVE
                B = complex(rand(n,d),rand(n,d));
                if rand(1)>0.5, B=real(B); end
                x1= tlsolve(G4,H4,B,[],xi,eta);
                T = toeplitz(c,r);
                x2= T\B;
                x3= tsolve(c,r,B);
                E = max([E,norm(x1-x2),norm(x1-x3)]);
            end
        end
    end
end

if abs(E)<value, s='ok'; else s='KO'; end
disp([s ' (' num2str(E) ')']);
warning(warn)

clear t t2 t3 t5 s s2 s3 s5 x1 x2 x3 xi xi1 eta eta1 G G1 G2 G3 G4 G5 H ...
    H1 H2 H3 H4 H5 Z_xi Z_eta c r B m n C T d E warn value
