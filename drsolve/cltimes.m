function y = cltimes(G,H,t,s,x)
%CLTIMES Cauchy-like matrix product.
%   Y = CLTIMES(G,H,t,s,X) computes the product C*X, where C is the
%   Cauchy-like matrix whose displacement equation is
%
%      diag(t) * C - C * diag(s) = G * H'
%
%   The matrix C must be fully-reconstructable, i.e., t(i)<>s(j) for
%   any (i,j), for the result to be correct.
%
%   See also ttimes, cl2full.

%   Antonio Arico' & Giuseppe Rodriguez, University of Cagliari, Italy
%   Email: {arico,rodriguez}@unica.it
%
%   Last revised Mar 25, 2010

C = intersect(t,s);
if numel(C)
    warning('drsolve:cltimes', 'C is not fully reconstructable');
end

m = size(G,1);
n = size(H,1);
t = t(:);
y = zeros(m,size(x,2));

for j=1:n
    %   x = x + (spdiags(t-s(j),0,m,m) \ (G*H(j,:)') )   * b(j,:); % NO
   y = y + ( (G*H(j,:)') ./ (t-s(j)) )   * x(j,:); % M=M+(v./w)*z 
end
