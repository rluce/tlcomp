function tests = toeppowersTest

tests = functiontests(localfunctions);

end

function test_zero(testCase)

c = 0;
r = 0;
s = 4;

[Gpowers, Bpowers] = toeppowers(c,r,s);
check_generators(testCase, zeros(1), Gpowers, Bpowers);

c = zeros(12,1);
r = zeros(1,12);
[Gpowers, Bpowers] = toeppowers(c,r,s);
check_generators(testCase, zeros(12), Gpowers, Bpowers);
end

function test_singleton(testCase)
c = 0.5;
r = 0.5;
s = 4;

T = 0.5;

[Gpowers, Bpowers] = toeppowers(c, r, s);
check_generators(testCase, T, Gpowers, Bpowers);

end

function test_jordan_block(testCase)
n = 9;

c = zeros(n,1);
r = zeros(1,n);
r(2) = 1.0;
s = 9;
T = toeplitz(c,r);
[Gpowers, Bpowers] = toeppowers(c, r, s);
check_generators(testCase, T, Gpowers, Bpowers);
end

function test_random(testCase)
n = 8;
s = 3;
[c,r] = random_toeplitz(n,n);

T = toeplitz(c,r);
nT = norm(T);
[Gpowers, Bpowers] = toeppowers(c, r, s);
check_generators(testCase, T, Gpowers, Bpowers, ...
    'AbsTol', nT*16*eps, 'RelTol', nT*16*eps);
end


function test_prolate(testCase)
% Use canonical generators
n = 64;
s = 13;
T = -gallery('prolate', n, 0.4);
c = T(:,1);
r = T(1,:);

[Gpowers, Bpowers] = toeppowers(c, r, s);
check_generators(testCase, T, Gpowers, Bpowers, 'AbsTol', 32*eps, 'RelTol', 32*eps);
end


function test_identity(testCase)

n = 13;
s = 20;

e1 = zeros(n,1);
e1(1) = 1;
c = e1;
r = e1;

T = toeplitz(c,r);
[Gpowers, Bpowers] = toeppowers(c, r, s);
check_generators(testCase, T, Gpowers, Bpowers);

end

function test_consistent_sizes(testCase)

[c,r] = lps_example3(16);
deg = 13;

[Gp, Bp] = toeppowers(c, r, deg);
for k=1:deg
    [~,s1] = size(Gp{k});
    [~,s2] = size(Bp{k});
    testCase.assertEqual(s1, s2);
    testCase.assertEqual(s1, 2*k);
end

end

function check_generators(testCase, T, Gpowers, Bpowers, varargin)
s = length(Bpowers);
T_true  = eye(size(T));

for i=1:s
    T_true = T_true * T;
    Tr = toeplkreconstruct(Gpowers{i}, Bpowers{i});
    testCase.assertEqual(Tr, T_true, varargin{:});
end


end
