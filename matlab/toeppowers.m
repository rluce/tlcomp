function [Gpowers, Bpowers] = toeppowers(c, r, s)
% [Bpowers, Tpowers] = power_generators(c, r, s)
%
% Compute generators for each monomial T^1, ..., T^s.

Gpowers = cell(1,s);
Bpowers = cell(1,s);

[G, B] = toepgen(c,r);

Gpowers{1} = G;
Bpowers{1} = B;

% Initilziation
% Caution: Here we use implicit knowledge about the structure of G and B.
Tg = G(:,1);
Tb = B(:,2);
Te1 = G(:,2);
Ten = B(:,1);

for i=2:s
    % At iteration i Tg = T^(i-1) * g, analoguous for the other three vectors
    tmp = toepmult(c, r, [Tg, Te1]);
    Tg  = tmp(:,1);
    Te1 = tmp(:,2);

    % At iteration i Tb = T'^(i-1) * b (note the ')
    tmp = toepmult(conj(r), conj(c), [Ten, Tb]);
    Ten = tmp(:,1);
    Tb  = tmp(:,2);
    
    % Extend G from power i-1 to power i
    G(:,end-1) = G(:,end-1) - Te1; % second last column needs update
    G = [G, Tg, Te1]; %#ok<AGROW>

    % Extend B from power i-1 to power i
    B(:,2) = B(:,2) - Ten;
    B = [Ten, Tb, B]; %#ok<AGROW>
    
    Gpowers{i} = G;
    Bpowers{i} = B;
end


end