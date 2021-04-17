classdef BitwiseSearch
    properties(Constant)
        Eps (1,1) double = 0.0001;
    end
        
    methods(Static)
        function [X, F] = Solve(a, b, f)
            arguments
                a (1,1) double
                b (1,1) double
                f handle
            end
         
            delta = (b - a) / 4.0;
            x0 = a; f0 = f(x0);
            x1 = x0; f1 = f0;
            
            while abs(delta) > Eps
                x0 = x1; f0 = f1;
                x1 = x0 + delta;
                f1 = f(x1);
                
                if (f1 < f0)
                    if (x1 > a && x1 < b) 
                        continue;
                    end
                end
                
                delta = -delta / 4.0;
            end
            
            X = x0; F = f0;
        end
    end
end

  
