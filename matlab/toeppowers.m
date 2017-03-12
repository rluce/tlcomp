function [Gpowers, Bpowers] = toeppowers(c, r, s, alg)
% [Bpowers, Tpowers] = power_generators(c, r, s, alg)
%
% Compute some generators for each power T^1, ... T^s.
%
% Paramter alg chooses between the 'full' and 'reduced' algorithm.

if nargin < 4
    alg = 'full';
end

switch(alg)
    case 'full'
        [Gpowers, Bpowers] = full_algorithm(c, r, s);
    case 'reduced'
        [Gpowers, Bpowers] = reduced_algorithm(c, r, s);
    otherwise
        error('expmt:InvalidParameter', ...
            'Invalid choice for alg parameter');
end

end

function [Gpowers, Bpowers] = reduced_algorithm(c,r,s)

if s < 1
    return;
end

Gpowers = cell(1,s);
Bpowers = cell(1,s);

[G, B] = stein_generator(c,r);

Gpowers{1} = G;
Bpowers{1} = B;

% Store intermediate factors
PG_times_G1 = G;
PB_times_B1 = B;
for i=2:s
    PG_times_G1 = apply_PG(c, r, PG_times_G1);
    PB_times_B1 = apply_PB(c, r, PB_times_B1);
    
    Gupdt = [ PG_times_G1, Gpowers{i-1}];
    % Implicit gauss elimination in order to make generators short.
    Gupdt(:,3) = Gupdt(:,3) - Gupdt(:,2);
    
    Bupdt = [ Bpowers{i-1}, PB_times_B1];
    Gpowers{i} = Gupdt;
    Bpowers{i} = Bupdt;
end


end

function [Gpowers, Bpowers] = full_algorithm(c,r,s)

if s < 1
    return;
end

Gpowers = cell(1,s);
Bpowers = cell(1,s);

[G, B] = stein_generator(c,r);

e1 = G(:,2);

Gpowers{1} = G;
Bpowers{1} = B;

% Store intermediate factors
PG_times_G1 = G;
PG_times_e1 = e1;
PB_times_B1 = B;
PB_times_e1 = e1;

for i=2:s
    PG_times_G1 = apply_PG(c, r, PG_times_G1);
    PG_times_e1 = apply_PG(c, r, PG_times_e1);

    PB_times_B1 = apply_PB(c, r, PB_times_B1);
    PB_times_e1 = apply_PB(c, r, PB_times_e1);

    Gupdt = [ PG_times_G1, Gpowers{i-1}, -PG_times_e1];
    Bupdt = [ Bpowers{i-1}(:,1:2*(i-1)), PB_times_B1, PB_times_e1, ...
        Bpowers{i-1}(:,end-i+3:end) ];
    Gpowers{i} = Gupdt;
    Bpowers{i} = Bupdt;
end

end

function Y = apply_PG(c, r, X)
Y = vapply(X, 'inv');
Y = toepmult(c,r,Y);
Y = vapply(Y);
end

function Y = apply_PB(c, r, X)
Y = vapply(X, 'inv');
Y = toepmult(conj(r), conj(c), Y);
Y = vapply(Y);
end

