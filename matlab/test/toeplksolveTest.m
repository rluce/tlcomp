classdef toeplksolveTest < matlab.unittest.TestCase
    properties (TestParameter)
        ctrans = {false, true};
    end
    
    methods (Test)
        
        function test_zerosolve(testCase, ctrans)
            G = 0;
            B = 0;
            
            b = 0;
            x_true = 0;
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, x_true);
            
            b = [0,0,0];
            x_true = [0,0,0];
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, x_true);
            
            G = [0,0,0,0,0];
            B = [0,0,0,0,0];
            b = [0,0,0];
            x_true = [0,0,0];
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, x_true);
        end
        
        function test_singleton(testCase, ctrans)
            
            G = 2;
            B = 1;
            
            b = 1;
            x_true = 1;
            x = toeplksolve(G, B, b, ctrans); % Same solution for A'
            testCase.assertEqual(x, x_true);
            
            b = [1,2,3];
            x_true = [1,2,3];
            x = toeplksolve(G, B, b, ctrans); % Same solution for A'
            testCase.assertEqual(x, x_true);
            
            G = 2i;
            B = 1;
            
            b = [1,2,3];
            
            if ~ctrans
                x_true = -1i * [1,2,3];
            else
                x_true = +1i * [1,2,3];
            end
            
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, x_true);
        end
        
        
        function test_identity(testCase, ctrans)
            n = 17;
            e1 = zeros(n,1);
            e1(1) = 1;
            [G, B] = toepgen(e1,e1);
            
            b = ones(n,1);
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, b, 'AbsTol', 4*eps, 'RelTol', 4*eps);

            b = (1:n)';
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, b, 'AbsTol', 16*eps, 'RelTol', 16*eps);
            
            b = randn(n,1);
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, b, 'AbsTol', 4*eps, 'RelTol', 4*eps);

            b = randn(n,1) + 1i * randn(n,1);
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, b, 'AbsTol', 4*eps, 'RelTol', 4*eps);

        end
        
        function test_prolate_matrix(testCase, ctrans)
            n = 17;
            T = gallery('prolate', n, 0.35);
            [G, B] = toepgen(T(:,1), T(1,:));
            b = [ones(n,1), linspace(0,1,n)'];
            x_true = T\b;
            x = toeplksolve(G,B,b, ctrans); % T is hermitian
            % Condition number around 10^5
            testCase.assertEqual(x, x_true, 'RelTol', 10^6 * n*eps);
            
        end
        
        function test_random_toeppoly(testCase, ctrans)
            n = 13;
            
            [c,r,T] = random_toeplitz(n,n);
            
            polycoef = [0.01, 0.1i, 1];
            pT = polyvalm(polycoef, T);
            [G, B] = toeppolyvalm(c,r,polycoef);
            
            b = [ones(n,1), linspace(0,1,n)'];
            if ~ctrans
                x_true = pT\b;
            else
                x_true = pT' \ b;
            end
            x = toeplksolve(G,B,b, ctrans);
            testCase.assertEqual(x, x_true, 'RelTol', 32*n*eps);
        end
    end
end
