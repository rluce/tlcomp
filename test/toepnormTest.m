classdef toepnormTest < matlab.unittest.TestCase
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
            val = toepnorm(c,c,p);
            testCase.assertEqual(val, 0);
            
            c = 1;
            val = toepnorm(c,c,p);
            testCase.assertEqual(val, 1);
            
            c = 1i;
            val = toepnorm(c,c,p);
            testCase.assertEqual(val, 1);
            
            c = rand;
            val = toepnorm(c,c,p);
            testCase.assertEqual(val, norm(c,p));
        end
        
        function test_zero(testCase, p)
            n = 16;
            c = zeros(n);
            val = toepnorm(c,c,p);
            testCase.assertEqual(val, 0);
        end
                
        function test_identity(testCase, p)
            n = 5;
            e1 = zeros(n,1);
            e1(1) = 1;
            val = toepnorm(e1,e1,p);
            testCase.assertEqual(val, norm(eye(n), p));
        end
        
        function test_random_real(testCase, p)
            n = 8;
            c = randn(n,1);
            r = randn(n,1);
            c(1) = r(1);
            val_true = norm(toeplitz(c,r), p);
            val = toepnorm(c,r,p);
            testCase.assertEqual(val, val_true, ...
                'RelTol', val_true * eps, 'AbsTol', val_true * eps);
        end

        function test_random_complex(testCase, p)
            [c,r,T] = random_toeplitz(9,9);
            val_true = norm(T, p);
            val = toepnorm(c,r,p);
            testCase.assertEqual(val, val_true, ...
                'RelTol', val_true * eps, 'AbsTol', val_true * eps);
        end

        
    end
end