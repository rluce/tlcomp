function tests = TLMatSteinTest

tests = functiontests(localfunctions);

end


function test_construct_empty(testCase)
TL = TLMatStein([]);
testCase.assertTrue(isempty(TL.G));
testCase.assertTrue(isempty(TL.B));

TL = TLMatStein([],[]);
testCase.assertTrue(isempty(TL.G));
testCase.assertTrue(isempty(TL.B));

TL = TLMatStein([],[], 'GB');
testCase.assertTrue(isempty(TL.G));
testCase.assertTrue(isempty(TL.B));

end


function test_construct_scalar(testCase)
A = 6;
TL = TLMatStein(A);
testCase.assertEqual(TL.G * TL.B', A);

G = 2;
B = 3;
TL = TLMatStein(G,B, 'GB');
testCase.assertEqual(TL.G * TL.B', A);

G = [1, 1];
B = [3, 3];
TL = TLMatStein(G,B, 'GB');
testCase.assertEqual(TL.G * TL.B', A);

c = 6;
r = 6;
TL = TLMatStein(c,r);
testCase.assertEqual(TL.G * TL.B', A);

c = 6;
r = 1;
testCase.assertError( @() TLMatStein(c,r), 'funmd:InconsistentInput');


end


function test_construct_eye(testCase)
n = 8;

E = eye(8);
e1 = E(:,1);
z = zeros(n,1);

TL = TLMatStein(E);
testCase.assertEqual(TL.G * TL.B', e1 * e1');

TL = TLMatStein(e1,e1);
testCase.assertEqual(TL.G * TL.B', e1 * e1');

TL = TLMatStein(e1,e1, 'GB');
testCase.assertEqual(TL.G * TL.B', e1 * e1');

TL = TLMatStein([e1,z], [e1,z], 'GB');
testCase.assertEqual(TL.G * TL.B', e1 * e1');

TL = TLMatStein([e1,z], [e1,z]);
testCase.assertEqual(TL.G * TL.B', e1 * e1');

end

function test_construct_random(testCase)
n = 9;
[c,r,T] = random_toeplitz(n,n);

TL = TLMatStein(T);
T2 = stein_reconstruction(TL.G, TL.B);
testCase.assertEqual(T2, T, 'RelTol', 100*eps);

TL = TLMatStein(c,r);
T2 = stein_reconstruction(TL.G, TL.B);
testCase.assertEqual(T2, T, 'RelTol', 100*eps);

[G, B] = stein_generator(c,r);
TL = TLMatStein(G,B);
T2 = stein_reconstruction(TL.G, TL.B);
testCase.assertEqual(T2, T, 'RelTol', 100*eps);

GG = [2 * G, -G];
BB = [B, B];
TL = TLMatStein(GG, BB);
T2 = stein_reconstruction(TL.G, TL.B);
testCase.assertEqual(T2, T, 'RelTol', 100*eps);


end

function test_construct_badinput(testCase)

testCase.assertError( @() TLMatStein(rand(4,2), rand(4,3)), ...
    'funmd:InconsistentInput');

end


function test_size(testCase)

TL = TLMatStein([]);
testCase.assertEqual(size(TL), [0,0]);

TL = TLMatStein(1);
testCase.assertEqual(size(TL), [1,1]);

TL = TLMatStein([1,2,5], [1,2,5]);
testCase.assertEqual(size(TL), [3,3]);

end