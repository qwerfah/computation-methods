function lab4()
    f = @(x) (tan((x.^4 + 2.*x.^2 - 2.*x + 2^0.5 + 1) ./ 8) + ...
          sin((4.*x.^3 - 7.*x - 9) ./ (20.*x + 28)));
      
    a = 0; b = 1;
    X = a:0.01:b;
    Y = f(X);
    eps = [0.01, 0.0001, 0.000001];
    
    NeutonMethodTable(f, a, b, X, Y, eps);
    FminbndTable(f, a, b, eps);
end

% Вывод таблицы для встроенной функции fminbnd
function FminbndTable(f, a, b, eps)
    fprintf("\nРезультаты вычисления точки минимума\n");
    fprintf("для различных значений точности функцией fminbnd:\n\n");
    fprintf(" # |    Eps   | N  |    x*   |   f(x*)\n");
    fprintf("---|----------|----|---------|--------\n");
    for i = 1:length(eps)
        options = optimset('TolX', eps(i));
        [x, fval, ~, output] = fminbnd(f, a, b, options);
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, eps(i), output.funcCount, x, fval);
    end
end

% Вывод таблицы и графиков для результатов работы метода Ньютона
function NeutonMethodTable(f, a, b, X, Y, eps)
    fprintf("\nРезультаты вычисления точки минимума\n");
    fprintf("для различных значений точности методом Ньютона:\n\n");
    fprintf("Eps - точность\n");
    fprintf("N - число обращений к целевой функции\n");
    fprintf("x* - найденная точка минимума функции\n");
    fprintf("f(x*) - найденный минимум функции\n\n");
    fprintf(" # |    Eps   | N  |    x*   |   f(x*)\n");
    fprintf("---|----------|----|---------|--------\n");
    
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    title('Метод Ньютона');
    
    for i = 1:length(eps)
        % Вычисление точки минимума и минимума функции
        [X0, F0, X_S, N] = NeutonMethod(a, b, f, eps(i));
        
        % Вывод строки таблицы результатов вычислений
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, eps(i), N, X0, F0);
        
        % Вывод графика для данной точности
        subplot(2, 2, i);
        % Целевая функция
        plot(X, Y, '-b','LineWidth',1.5);
        hold on;
        % Последовательность приближений
        plot(X_S, f(X_S), '*g','LineWidth', 2);
        % Точка минимума
        plot(X0, F0, '*r','LineWidth', 4);
        title(sprintf("Точность Eps = %2.0e", eps(i)));
        legend('Целевая функция', 'Последовательность приближений', ...
            'Точка минимума');
    end
end

% Разностная аппроксимация первой производной
function df = df(f_prev, f_next, delta)
    df = (f_next - f_prev) / 2 / delta;
end

% Разностная аппроксимация второй производной
function d2f = d2f(f_prev, f_m, f_next, delta)
    d2f = (f_next - 2*f_m + f_prev) / delta / delta;
end

% Метод Ньютона
function [X, F, X0, N] =  NeutonMethod(a, b, f, eps)
    x0_prev = (a + b) / 1.3;
    x_prev = x0_prev - eps;
    x_next = x0_prev + eps;
    
    f_prev = f(x_prev);
    f_next = f(x_next);
    
    d2ff = d2f(f_prev, f(x0_prev), f_next, x_next - x_prev);

    X0 = [x0_prev];
    N = 3; % Число обращений к целевой функции
    
    while true
        % Очередное приближение точки минимума
        dff = df(f_prev, f_next, x_next - x_prev);
        x0 = x0_prev - dff / d2ff;
        
        % Проверка условия завершения поиска
        if abs(x0_prev - x0) <= eps
            break;
        end
        
        % Вычисление точек для разностной аппроксимации производной
        x0_prev = x0;
        x_prev = x0 - eps; 
        x_next = x0 + eps;
        f_prev = f(x_prev);
        f_next = f(x_next);
        N = N + 2;
        X0 = [X0, x0];
    end
    
    X = x0;
    F = f_next;
end