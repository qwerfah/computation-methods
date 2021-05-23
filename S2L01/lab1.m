function lab1()
    f = @(x) (tan((x.^4.0 + 2.0.*x.^2.0 - 2.0.*x + 2.0^(1.0/2.0) + 1.0) ./ 8.0) + ...
          sin((4.0.*x.^3.0 - 7.0.*x - 9.0) ./ (20.0.*x + 28.0)));
      
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
        [X0, F0, X_S, F_S, N] = BitwiseSearch(0, 1, f, eps(i));
        
        % Вывод строки таблицы результатов вычислений
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, eps(i), N, X0, F0);
        
        % Вывод графика для данной точности
        ax = nexttile;
        % Целевая функция
        plot(ax, X, Y, '-b','LineWidth',1.5);
        hold on;
        % Последовательность приближений
        plot(ax, X_S, F_S, '*g','LineWidth', 2);
        % Точка минимума
        plot(ax, X0, F0, '*r','LineWidth', 4);
        title(ax, sprintf("Точность Eps = %2.0e", eps(i)));
        legend('Целевая функция', 'Последовательность приближений', ...
            'Точка минимума');
    end
    
end

% Метод поразрядного поиска
function [X, F, X_S, F_S, N] = BitwiseSearch(a, b, f, eps)
    arguments
        a   double           % Левая граница отрезка
        b   double           % Левая граница отрезка
        f   function_handle  % Целевая функция
        eps double           % Точность
    end
    
    delta = (b - a) / 4.0;  % Шаг поиска

    x0 = a;      % Начальное приближение точки минимума
    f0 = f(x0);  % Начальное приближение минимума функции
    N = 1;
    
    X_S = [];    % Массив всех приближений точки минимума
    F_S = [];    % Массив всех приближений минимума функции

    x1 = x0; 
    f1 = f0;

    while abs(delta) > eps
        X_S = [X_S, x1];
        F_S = [F_S, f1];

        x0 = x1;
        f0 = f1;

        x1 = x0 + delta;
        f1 = f(x1);
        N = N + 1;

        if (f1 < f0 && x1 > a && x1 < b)
            continue;
        end

        delta = -delta / 4.0;
    end

    X = x1; 
    F = f1;
end

  
