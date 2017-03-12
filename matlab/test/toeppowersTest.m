function tests = toeppowersTest

tests = functiontests(localfunctions);

end

function test_singleton(testCase)
c = 0.5;
r = 0.5;
s = 4;

try_all_algs(testCase, c, r, s);

end

function test_reduced_size(testCase)
% Check generator sizes of reduced algorithm
n = 17;
s = 5;
[c,r] = random_toeplitz(n,n);
[Gpowers, Bpowers] = toeppowers(c, r, s, 'reduced');

for i=1:s
    testCase.assertEqual(size(Gpowers{i}), [n,2*i]);
    testCase.assertEqual(size(Bpowers{i}), [n,2*i]);
end

end

function test_random(testCase)
% Use canonical generators
n=8;
s = 3;
[c,r] = random_toeplitz(n,n);
try_all_algs(testCase, c, r, s);
end


function test_prolate(testCase)
% Use canonical generators
n = 64;
s = 13;
T = -gallery('prolate', n, 0.4);
try_all_algs(testCase, T(:,1), T(1,:), s);
end


function test_identity(testCase)

n = 13;
s = 20;

e1 = zeros(n,1);
e1(1) = 1;
try_all_algs(testCase, e1, e1, s);


end

function test_consistent_sizes(testCase)

[c,r] = lps_example3(16);
deg = 13;

for alg = {'full', 'reduced'}
    [Gp, Bp] = toeppowers(c, r, deg, alg{:});
    for k=1:13
        [~,s1] = size(Gp{k});
        [~,s2] = size(Bp{k});
        testCase.assertEqual(s1, s2);
    end

end

end

function try_all_algs(testCase, c, r, s)

T = toeplitz(c,r);

[Gpowers, Bpowers] = toeppowers(c, r, s, 'full');
check_generators(testCase, T, Gpowers, Bpowers);

[Gpowers, Bpowers] = toeppowers(c, r, s, 'reduced');
check_generators(testCase, T, Gpowers, Bpowers);


end

function check_generators(testCase, T, Gpowers, Bpowers)
n = size(T,1);
s = length(Bpowers);
T_true  = eye(size(T));
for i=1:s
    T_true = T_true * T;
    Tr = stein_reconstruction(Gpowers{i}, Bpowers{i});
    testCase.assertEqual(Tr, T_true, 'RelTol', 64*n*eps, 'AbsTol', 64*n*eps);
end


end
