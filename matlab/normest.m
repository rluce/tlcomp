function [e,cnt] = normest(S, tol)
%NORMEST Estimate the matrix 2-norm.


if nargin < 2
    tol = 1.e-6;
end
maxiter = 100;

x = ones(size(S,1), 1);
cnt = 0;
e = norm(x);
if e == 0
    return;
end
x = x/e;
e0 = 0;
while abs(e-e0) > tol*e
   e0 = e;
   Sx = S*x;
   if nnz(Sx) == 0
      Sx = rand(size(Sx),class(Sx));
   end
   x = S'*Sx;
   normx = norm(x);
   e = normx/norm(Sx);
   x = x/normx;
   cnt = cnt+1;
   if cnt > maxiter
      warning(message('tlcomp:normest:notconverge', maxiter, sprintf('%g',tol)));
      break;
   end
end
