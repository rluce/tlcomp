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
            testCase.assertEqual(x, b, 'AbsTol', 32*eps, 'RelTol', 16*eps);

            b = (1:n)';
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, b, 'AbsTol', 128*eps, 'RelTol', 128*eps);
            
            b = randn(n,1);
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, b, 'AbsTol', 32*eps, 'RelTol', 16*eps);

            b = randn(n,1) + 1i * randn(n,1);
            nb = norm(b);
            x = toeplksolve(G, B, b, ctrans);
            testCase.assertEqual(x, b, 'AbsTol', 8*nb*eps, 'RelTol', 8*nb*eps);

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
            nfact = cond(pT) * norm(b);
            testCase.assertEqual(x, x_true, ...
                'AbsTol', 4*nfact*eps, 'RelTol', 4*nfact*eps);
        end
        
        function test_default_notranspose(testCase)
            n = 7;
            c = 1:n;
            r = zeros(n,1);
            r(1) = 1;
            
            T = toeplitz(c,r);
            [G, B] = toepgen(c,r);
            b = ones(n,1);
            x_true = zeros(n,1);
            x_true(1) = 1;
            x_true(2) = -1;   
            x = toeplksolve(G,B,b); % No value for ctrans, so 'no transpose'
            testCase.assertEqual(x, x_true, 'RelTol', 4*n*eps, 'AbsTol', 4*n*eps);

           
        end
    end
end
