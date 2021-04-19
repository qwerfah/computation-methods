classdef GoldenRatio
    properties(Constant)
        Eps (1,1) double = 0.000001;
        Tau (1,1) double = 0.61803;
    end
    
    methods(Static)
        function [X, F, A_S, B_S] = Solve(a, b, f)
            arguments
                a (1,1) double
                b (1,1) double
                f function_handle
            end
            
            A_S = [];
            B_S = [];
            
            l = b - a;
            
            x1 = b - GoldenRatio.Tau * l; 
            x2 = a + GoldenRatio.Tau * l; 
            
            f1 = f(x1);
            f2 = f(x2);
            
            while true
                A_S = [A_S, a];
                B_S = [B_S, b];
                
                if (f1 >= f2)
                    a = x1;
                    l = b - a;
                    
                    x1 = x2; 
                    f1 = f2;
                    
                    x2 = a + GoldenRatio.Tau * l;
                    f2 = f(x2);
                else
                    b = x2;
                    l = b - a;
                    
                    x2 = x1;
                    f2 = f1;
                    
                    x1 = b - GoldenRatio.Tau * l;
                    f1 = f(x1);
                end
                
                if (l < 2 * GoldenRatio.Eps)
                    break; 
                end
            end
            
            X = (a + b) / 2.0;
            F = f(X);
        end
    end
end