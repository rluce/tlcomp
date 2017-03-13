function tests = fcirculantTest

tests = functiontests(localfunctions);

end

function test_singleton(testCase)

n = 1;
f = 0;
Z = fcirculant(n,f);
testCase.assertEqual(Z, sparse(1,1));

n = 1;
f = 1;
Z = fcirculant(n,f);
testCase.assertEqual(Z, spones(1));

n = 1;
f = -1;
Z = fcirculant(n,f);
testCase.assertEqual(Z, -spones(1));


end

function test_fzero(testCase)

n = 17;
Z = fcirculant(n, 0);
Z_true = downshift(n);
testCase.assertEqual(Z, Z_true);

end

function test_fnonzero(testCase)

n = 9;
f = pi;
Z = fcirculant(n, f);
ZZ = downshift(n);
ZZZ = sparse(n,n);
ZZZ(1,n) = f;

testCase.assertEqual(Z, ZZ+ZZZ);

end

