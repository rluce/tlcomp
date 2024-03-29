classdef toeplkmultTest < matlab.unittest.TestCase
    properties (TestParameter)
        alg = {'full', 'fft_naive', 'fft'};
        ctrans = {false, true};
    end

    methods (Test)
        function test_singleton(testCase, ctrans, alg)
            c = 0;
            r = 0;
            [G,B] = toepgen(c,r);
            
            x = 1;
            y_true = 0;
            y = toeplkmult(G, B, x, ctrans, alg);
            testCase.assertEqual(y, y_true, 'AbsTol', eps);
            
            c = 1;
            r = 1;
            [G,B] = toepgen(c,r);
            x = 1;
            y_true = 1;
            y = toeplkmult(G, B, x, ctrans, alg);
            testCase.assertEqual(y, y_true, 'AbsTol', eps);
            
            c = 1i;
            r = 1i;
            [G,B] = toepgen(c,r);
            
            x = 1i;
            
            if ~ctrans
                y_true = -1;
            else
                y_true = 1;
            end
            
            y = toeplkmult(G, B,  x, ctrans, alg);
            testCase.assertEqual(y, y_true, 'AbsTol', eps);
        end

        function test_identity(testCase, ctrans, alg)
            n = 13;
            
            e1 = zeros(n,1);
            e1(1) = 1.0;
            
            c = e1;
            r = e1;
            [G,B] = toepgen(c,r);
            
            x = e1;
            y_true = e1;
            y = toeplkmult(G,B,x,ctrans,alg);
            testCase.assertEqual(y, y_true, 'AbsTol', n*eps);
            
            x = randn(n,1);
            y_true = x;
            y = toeplkmult(G,B,x,ctrans,alg);
            testCase.assertEqual(y, y_true, 'AbsTol', n*eps);
            
            x = randn(n,1) + 1i * randn(n,1);
            y_true = x;
            y = toeplkmult(G,B,x,ctrans,alg);
            testCase.assertEqual(y, y_true, 'AbsTol', n*eps);
        end
        
        function test_tril_ones(testCase, ctrans, alg)
            n = 15;
            
            c = ones(n,1);
            r = zeros(n,1);
            r(1) = 1;
            [G,B] = toepgen(c,r);
            
            x = ones(n,1);
            
            if ~ctrans
                y_true = (1:n)';
            else
                y_true = (n:-1:1)';
            end
            
            y = toeplkmult(G,B,x,ctrans,alg);
            testCase.assertEqual(y, y_true, 'RelTol', 4*n*eps);
        end
        
        function test_triltoep_real(testCase, ctrans, alg)
            n = 14;
            c = randn(n,1);
            r = zeros(n,1);
            r(1) = c(1);
            [G,B] = toepgen(c,r);
            nT = norm(toeplitz(c,r));
            
            x = ones(n,1);
            y_true = cumsum(c);
            if ctrans
                y_true = y_true(end:-1:1);
            else
                
            end
            y = toeplkmult(G,B,x,ctrans,alg);
            testCase.assertEqual(y, y_true, 'AbsTol', 16*nT*eps, 'RelTol', 16*nT*eps);
            
        end

        function test_triltoep_complex(testCase, ctrans, alg)
            n = 14;
            c = randn(n,1) + 1i * randn(n,1);
            r = zeros(n,1);
            r(1) = c(1);
            [G,B] = toepgen(c,r);
            
            x = ones(n,1);
            y_true = cumsum(c);
            if ctrans
                y_true = conj(y_true(end:-1:1));
            else
                
            end
            y = toeplkmult(G,B,x,ctrans,alg);
            testCase.assertEqual(y, y_true, 'RelTol', 256*n*eps);
            
        end

        
        function test_random_real_rank1(testCase, ctrans, alg)
            n = 9;
            G = orth(randn(n,1));
            B = orth(randn(n,1));
            T = toeplkreconstruct(G,B);
            x = linspace(1,2,n)';
            if ~ctrans
                y_true = T * x;
            else
                y_true = T' * x;
            end
            y = toeplkmult(G,B,x,ctrans,alg);
            testCase.assertEqual(y, y_true, 'RelTol', n*eps, 'AbsTol', n*eps);
            
            testCase.assertTrue(isreal(y));
        end
        
        
        function test_random_real_rank4(testCase, ctrans, alg)
            n = 9;
            G = orth(randn(n,4));
            B = orth(randn(n,4));
            T = toeplkreconstruct(G,B);
            x = ones(n,1);
            
            if ~ctrans
                y_true = T*x;
            else
                y_true = T'*x;
            end
            y = toeplkmult(G,B,x,ctrans,alg);
            testCase.assertEqual(y, y_true, 'AbsTol', 128*n*eps, 'RelTol', 128*n*eps);
            testCase.assertTrue(isreal(y));

        end
        
        function test_random_complex_rank4(testCase, ctrans, alg)
            n = 9;
            G = orth(randn(n,4) + 1i * randn(n,4));
            B = orth(randn(n,4) + 1i * randn(n,4));
            T = toeplkreconstruct(G,B);
            x = ones(n,1);
            if ~ctrans
                y_true = sum(T, 2);
            else
                y_true = sum(T, 1)';
            end
            y = toeplkmult(G,B,x,ctrans,alg);
            nfact = 4*norm(T);
            testCase.assertEqual(y, y_true, 'RelTol', nfact*eps, 'AbsTol', nfact*eps);
        end
        
        
        function test_invtoep(testCase, ctrans, alg)
            n = 14;
            c = rand(n,1);
            r = rand(n,1);
            c(1) = c(1) + n;
            r(1) = c(1);
            T = toeplitz(c,r);
            e = ones(n,1);
            if ~ctrans
                x = T*e;
            else
                x = T' * e;
            end
            
            [Ginv, Binv] = toepinvgen(c,r);
            y = toeplkmult(Ginv, Binv, x, ctrans, alg);
            testCase.assertEqual(y, e, 'AbsTol', 8192*n*eps);
        end
        
        
        function test_squared_toep(testCase, ctrans, alg)
            n=8;
            T = gallery('prolate', n, 0.48);
            c = T(:,1);
            r = T(1,:)';
            [G,B] = toepgen(c,r);
            [Gs, Bs] = toeplksquare(G, B);
            x = ones(n,1);
            if ~ctrans
                y_true = T*(T*x);
            else
                y_true = (T*T)' * x;
            end
            y = toeplkmult(Gs, Bs, x, ctrans, alg);
            testCase.assertEqual(y, y_true, 'RelTol', n*eps);
        end
        
        
        function test_sum_toep(testCase, ctrans, alg)
            n = 15;
            c1 = orth(randn(n,1));
            r1 = orth(randn(n,1));
            r1(1) = c1(1);
            c2 = orth(randn(n,1));
            r2 = orth(randn(n,1));
            r2(1) = c2(1);
            
            [G1, B1] = toepgen(c1,r1);
            [G2, B2] = toepgen(c2,r2);
            G = [G1, G2];
            B = [B1, B2];
            T = toeplitz(c1+c2, r1+r2);
            
            x = rand(n,1);
            if ~ctrans
                y_true = T*x;
            else
                y_true = T' * x;
            end
            y = toeplkmult(G, B, x, ctrans, alg);
            testCase.assertEqual(y, y_true, 'AbsTol', n*eps);
        end
        
        function test_multrhs(testCase, ctrans, alg)
            
            n = 7;
            
            G = orth(randn(n,3));
            B = orth(randn(n,3));
            T = toeplkreconstruct(G,B);
            
            x = eye(n,2);
            y = toeplkmult(G,B,x, ctrans, alg);
            nfact = 4 * norm(T) * norm(y);
            if ~ctrans
                testCase.assertEqual(y(:,1), T(:,1), ...
                    'RelTol', nfact * eps, 'AbsTol', nfact * eps);
                testCase.assertEqual(y(:,2), T(:,2), ...
                    'RelTol', nfact * eps, 'AbsTol', nfact * eps);
            else
                testCase.assertEqual(y(:,1), T(1,:)', ...
                    'RelTol', nfact * eps, 'AbsTol', nfact * eps);
                
                testCase.assertEqual(y(:,2), T(2,:)', ...
                    'RelTol', nfact * eps, 'AbsTol', nfact * eps);
            end
        end
    end
end
