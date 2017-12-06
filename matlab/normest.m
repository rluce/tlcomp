function [e,cnt] = normest(TL, tol)
%NORMEST Estimate the matrix 2-norm.


if nargin < 2
    tol = 1.e-6;
end
maxiter = 100;

x = ones(size(TL,1), 1);
cnt = 0;
e = norm(x);
x = x/e;
e0 = 0;
while abs(e-e0) > tol*e
   e0 = e;
   Sx = TL*x;
   if nnz(Sx) == 0
      Sx = rand(size(Sx),class(Sx));
   end
   x = TL'*Sx;
   normx = norm(x);
   e = normx/norm(Sx);
   x = x/normx;
   cnt = cnt+1;
   if cnt > maxiter
      warning('tlcomp:normest:notconverge', 'normest did not converge');
      break;
   end
end
