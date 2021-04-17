method = BitwiseSearch();

method.Solve(-5.0, 5.0, @f);

function F = f(x)
    arguments
        x (1, 1) double
    end
    
    F = x * x;
end
