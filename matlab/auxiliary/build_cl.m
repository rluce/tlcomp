function C = build_cl(s, t, G, B)
% C = build_cl(s, t, G, B) -- construct CL matrix from data

s = s(:);
t = t(:);

m = length(s);
n = length(t);

C = G*B ./ ([s, -ones(m,1)] * [ones(1,n); transpose(t)]);

end % of function
