classdef toeplkprodTest < matlab.unittest.TestCase
    properties (TestParameter)
        alg = {'full', 'fft'};
    end
    
    methods (Test)
        
        function test_eye(testCase, alg)
            G1 = 2;
            B1 = 1;
            G2 = 2;
            B2 = 1;
            
            [GP, BP] = toeplkprod(G1, B1, G2, B2, alg);
            testCase.assertEqual(toeplkreconstruct(GP, BP), 1);
            
            e1 = zeros(8,1);
            e1(1) = 1;
            e8 = zeros(8,1);
            e8(8) = 1;
            G1 = 2*e1;
            B1 = e8;
            G2 = 2*e1;
            B2 = e8;
            
            [GP, BP] = toeplkprod(G1, B1, G2, B2, alg);
            testCase.assertEqual(toeplkreconstruct(GP, BP), eye(8));
        end
        
        function test_square(testCase, alg)
            n = 12;
            G1 = randn(n,4);
            B1 = randn(n,4);
            G2 = G1;
            B2 = B1;
            P_true = toeplkreconstruct(G1, B1)^2;
            [GP, BP] = toeplkprod(G1, B1, G2, B2, alg);
            P = toeplkreconstruct(GP, BP);
            maxP = max(abs(P_true(:)));
            testCase.assertEqual(P, P_true, 'AbsTol', 16*maxP*n*eps, 'RelTol', 16*n*eps);
        end
        
        function test_random(testCase, alg)
            G1 = randn(12,4);
            B1 = randn(12,4);
            G2 = 1i * randn(12,3);
            B2 = 1i * randn(12,3);
            P_true = toeplkreconstruct(G1, B1) * toeplkreconstruct(G2, B2);
            [GP, BP] = toeplkprod(G1, B1, G2, B2, alg);
            P = toeplkreconstruct(GP, BP);
            maxP = max(abs(P_true(:)));
            testCase.assertEqual(P, P_true, 'AbsTol', maxP*16*eps, 'RelTol', 16*eps);
        end
        
        function test_toeplitz(testCase, alg)
            [c1, r1, T1] = random_toeplitz(8,8);
            [c2, r2, T2] = random_toeplitz(8,8);
            [G1, B1] = toepgen(c1,r1);
            [G2, B2] = toepgen(c2,r2);
            [GP, BP] = toeplkprod(G1, B1, G2, B2, alg);
            P = toeplkreconstruct(GP, BP);
            P_true = T1 * T2;
            maxP = max(abs(P_true(:)));
            testCase.assertEqual(P, T1 * T2, 'AbsTol', 16*maxP*eps, 'RelTol', 16*eps);
        end
    end
end

