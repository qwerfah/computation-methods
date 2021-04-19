classdef BitwiseSearch
    properties(Constant)
        Eps (1,1) double = 0.000001;
    end
        
    methods(Static)
        function [X, F, X_S] = Solve(a, b, f)
            arguments
                a (1,1) double
                b (1,1) double
                f function_handle
            end
            
            X_S = [];
            
            delta = (b - a) / 4.0;
            
            x0 = a; 
            f0 = f(x0);
            
            x1 = x0; 
            f1 = f0;
            
            while abs(delta) > BitwiseSearch.Eps
                X_S = [X_S, x1];
                
                x0 = x1; 
                f0 = f1;
                
                x1 = x0 + delta;
                f1 = f(x1);
                
                if (f1 < f0 && x1 > a && x1 < b)
                    continue;
                end
                
                delta = -delta / 4.0;
            end
            
            X = x1; 
            F = f1;
        end
    end
end

  
