function lab3()
    f = @(x) (tan((x.^4 + 2.*x.^2 - 2.*x + 2^0.5 + 1) ./ 8) + ...
          sin((4.*x.^3 - 7.*x - 9) ./ (20.*x + 28)));
      
    a = 0; b = 1;
    X = a:0.01:b;
    Y = f(X);
    eps = [0.01, 0.0001, 0.000001];
    
    fprintf("\nРезультаты вычисления точки минимума\n");
    fprintf("для различных значений точности:\n\n");
    fprintf("Eps - точность\n");
    fprintf("N - число обращений к целевой функции\n");
    fprintf("x* - найденная точка минимума функции\n");
    fprintf("f(x*) - найденный минимум функции\n\n");
    fprintf(" # |    Eps   | N  |    x*   |   f(x*)\n");
    fprintf("---|----------|----|---------|--------\n");
    
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    title('Метод поразрядного поиска');
    tiledlayout(2, 2);
    
    for i = 1:length(eps)
        % Вычисление точки минимума и минимума функции
        [X0, F0, X_S, N] = NeutonMethod(a, b, f, eps(i));
        
        % Вывод строки таблицы результатов вычислений
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, eps(i), N, X0, F0);
        
        % Вывод графика для данной точности
        ax = nexttile;
        % Целевая функция
        plot(ax, X, Y, '-b','LineWidth',1.5);
        hold on;
        % Последовательность приближений
        plot(ax, X_S, f(X_S), '*g','LineWidth', 2);
        % Точка минимума
        plot(ax, X0, F0, '*r','LineWidth', 4);
        title(ax, sprintf("Точность Eps = %2.0e", eps(i)));
        legend('Целевая функция', 'Последовательность приближений', ...
            'Точка минимума');
    end
    
end

function df = df(f_prev, f_next, delta)
    df = (f_next - f_prev) / 2 / delta;
end

function d2f = d2f(f_prev, f_m, f_next, delta)
    d2f = (f_next - 2*f_m + f_prev) / delta / delta;
end

function [X, F, X0, N] =  NeutonMethod(a, b, f, eps)
    arguments
        a    double            % Левая граница отрезка
        b    double            % Правая граница отрезка
        f    function_handle   % Целевая функция
        eps  double            % Точность
    end
    
    x0 = a; x_prev = a; x_next = b;
    f_prev = f(x_prev); f_next = f(x_next); f_m = f(x0);
    
    % f1 = f(a); f2 = f(b); 
    X0 = [x0];
    n = 3; % Число обращений к целевой функции
    
    while true
        dff = df(f_prev, f_next, abs(x_next - x_prev));
        d2ff = d2f(f_prev, f_m, f_next, abs(x_next - x_prev));
        x0 = x0 - dff / d2ff;
        
        f_prev = f_m;
        f_m = f_next;
        f_next = f(x0);
        n = n + 1;
        
        if abs(df(f_prev, f_next, abs(x_next - x_prev))) <= eps
            break;
        end
        
        
        X0 = [X0, x0];
        
    end
    
    X = x0;
    F = f_next;
    N = n;
end