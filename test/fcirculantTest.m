function tests = fcirculantTest

tests = functiontests(localfunctions);

end

function test_empty(testCase)

testCase.assertTrue(isempty(fcirculant(0)));
testCase.assertTrue(isempty(fcirculant(0,1)));
testCase.assertTrue(isempty(fcirculant(0,-2i)));

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
Z_true = spdiags(ones(n,1), -1, n, n);
testCase.assertEqual(Z, Z_true);

end

function test_fnonzero(testCase)

n = 9;
f = pi;
Z = fcirculant(n, f);
ZZ = spdiags(ones(n,1), -1, n, n);
ZZZ = sparse(n,n);
ZZZ(1,n) = f;

testCase.assertEqual(Z, ZZ+ZZZ);

end


function test_issparse(testCase)
% In several places we need to apply the unit f-ciruclant matrix Z to a
% vector or a matrix.  This amounts to scaling and permuting the operand,
% but we are lazy to do that on our own.  If Z is in sparse format,
% Matlab's sparse matrix multiply achieves already the best complexity for
% this permutation, so we simply resort to Z*A.  Hence, we really need
% that Z is in sparse format.

n = 83;
Z0 = fcirculant(n);
Z1 = fcirculant(n, 1);
Zm1 = fcirculant(n, -1);
Zpi = fcirculant(n, pi);

testCase.assertTrue(issparse(Z0));
testCase.assertTrue(issparse(Z1));
testCase.assertTrue(issparse(Zm1));
testCase.assertTrue(issparse(Zpi));
end
