classdef toeplksquareTest < matlab.unittest.TestCase
    properties (TestParameter)
        alg = {'full', 'fft'};
    end
    
    methods (Test)
        
        function test_singleton(testCase, alg)
            % 1x1 zero
            G = 0;
            B = 0;
            [Gs, Bs] = toeplksquare(G, B);
            T2 = toeplkreconstruct(Gs, Bs);
            testCase.assertEqual(T2, 0);
            
            % 1x1 identity
            G = 2;
            B = 1;
            Ts_true = 1.0;
            
            [Gs, Bs] = toeplksquare(G, B, alg);
            Ts = toeplkreconstruct(Gs, Bs);
            testCase.assertEqual(Ts, Ts_true, 'AbsTol', eps);
        end
        
        function test_identity(testCase, alg)
            n=15;
            e1 = zeros(n, 1);
            e1(1) = 1.0;
            [G, B] = toepgen(e1, e1);
            Ts_true = eye(n);

            [Gs, Bs] = toeplksquare(G, B, alg);
            Ts = toeplkreconstruct(Gs, Bs);
            testCase.assertEqual(Ts, Ts_true, 'AbsTol', n*eps);
        end
        
        function test_random_real(testCase, alg)
            n = 12;
            c = randn(n,1);
            r = randn(n,1);
            r(1) = c(1);
            T = toeplitz(c,r);
            Ts_true = T*T;
            [G, B] = toepgen(c,r);

            [Gs, Bs] = toeplksquare(G, B, alg);
            Ts = toeplkreconstruct(Gs, Bs);
            testCase.assertEqual(Ts, Ts_true, 'AbsTol', 4*n*eps, 'RelTol', 4*n*eps);
        end
        
        function test_random_complex(testCase, alg)
            n = 12;
            c = randn(n,1) + 1i*randn(n,1);
            r = randn(n,1) + 1i*randn(n,1);
            r(1) = c(1);
            T = toeplitz(c,r);
            Ts_true = T*T;
            [G, B] = toepgen(c,r);
            [Gs, Bs] = toeplksquare(G, B, alg);
            Ts = toeplkreconstruct(Gs, Bs);
            testCase.assertEqual(Ts, Ts_true, 'AbsTol', 4*n*eps, 'RelTol', 4*n*eps);
        end
        
        function test_random_d5(testCase, alg)
            n = 16;
            d = 5;
            G = randn(n,d) + 1i * randn(n,d);
            B = randn(n,d) + 1i * randn(n,d);
            T = toeplkreconstruct(G, B);
            Ts_true = T * T;
            [Gs, Bs] = toeplksquare(G, B, alg);
            Ts = toeplkreconstruct(Gs, Bs);
            testCase.assertEqual(Ts, Ts_true, 'RelTol', 64*n*eps, 'AbsTol', 64*n*eps);
            
        end
    end
end


