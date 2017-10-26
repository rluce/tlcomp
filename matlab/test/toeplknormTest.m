classdef toeplknormTest < matlab.unittest.TestCase
    properties (TestParameter)
        p = {1, inf, 'fro'};
    end

    methods (Test)
        
        function test_empty(testCase, p)
            val = toeplknorm([], [], p);
            testCase.assertEqual(val, 0);
        end
        
        function test_singleton(testCase, p)
            c = 0;
            [G, B] = toepgen(c,c);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, 0);
            
            c = 1;
            [G, B] = toepgen(c,c);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, 1);
            
            c = 1i;
            [G, B] = toepgen(c,c);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, 1);
            
            c = rand;
            [G, B] = toepgen(c,c);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, norm(c,p));
        end
        
        function test_zero(testCase, p)
            n = 16;
            c = zeros(n);
            [G, B] = toepgen(c,c);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, 0);
            
            G = zeros(n,5);
            B = zeros(n,5);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, 0);
        end
                
        function test_identity(testCase, p)
            n = 5;
            e1 = zeros(n,1);
            e1(1) = 1;
            [G, B] = toepgen(e1,e1);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, norm(eye(n), p));
        end
        
        function test_random_toep_real(testCase, p)
            n = 8;
            c = randn(n,1);
            r = randn(n,1);
            c(1) = r(1);
            val_true = norm(toeplitz(c,r), p);
            [G, B] = toepgen(c,r);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, val_true, ...
                'RelTol', val_true * eps, 'AbsTol', val_true * eps);
        end

        function test_random_toeplk_real(testCase, p)
            n = 12;
            d = 4;
            G = randn(n,d);
            B = randn(n,d);
            val_true = norm(toeplkreconstruct(G, B), p);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, val_true, ...
                'RelTol', val_true * eps, 'AbsTol', val_true * eps);
        end

        
        function test_random_toeplk_complex(testCase, p)
            n = 11;
            d = 5;
            G = randn(n,d) + 1i * randn(n,d);
            B = randn(n,d) + 1i * randn(n,d);
            val_true = norm(toeplkreconstruct(G, B), p);
            val = toepnorm(G,B,p);
            testCase.assertEqual(val, val_true, ...
                'RelTol', val_true * eps, 'AbsTol', val_true * eps);
        end

        
    end
end