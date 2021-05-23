function lab3()
    f = @(x) (tan((x.^4 + 2.*x.^2 - 2.*x + 2^0.5 + 1) ./ 8) + ...
          sin((4.*x.^3 - 7.*x - 9) ./ (20.*x + 28)));
      
    a = 0; b = 1;
    X = a:0.01:b;
    Y = f(X);
    eps = [0.01, 0.0001, 0.000001];
    
    NeutonMethodTable(f, a, b, X, Y, eps);
end

% Вывод сравнительной таблицы методов для точности 10^-6
function ComparisonTable(f, a, b, X, Y)
    arguments
        f   function_handle  % Целевая функция
        a   double           % Левая граница отрезка
        b   double           % Правая граница отрезка
        X   (1,:) double     % Массив значений аргумента целевой функции
        Y   (1,:) double     % Массив значений целевой ыункции
    end
    
    fprintf(" # |      Метод     | N  |    x*   |   f(x*)\n");
    fprintf("---|----------------|----|---------|--------\n");
end

% Вывод таблицы и графиков для результатов работы метода Ньютона
function NeutonMethodTable(f, a, b, X, Y, eps)
    arguments
        f   function_handle  % Целевая функция
        a   double           % Левая граница отрезка
        b   double           % Правая граница отрезка
        X   (1,:) double     % Массив значений аргумента целевой функции
        Y   (1,:) double     % Массив значений целевой ыункции
        eps double           % Массив значений точности
    end

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
    
    for e = eps
        % Вычисление точки минимума и минимума функции
        [X0, F0, X_S, N] = NeutonMethod(a, b, f, e);
        
        % Вывод строки таблицы результатов вычислений
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, e, N, X0, F0);
        
        % Вывод графика для данной точности
        ax = nexttile;
        % Целевая функция
        plot(ax, X, Y, '-b','LineWidth',1.5);
        hold on;
        % Последовательность приближений
        plot(ax, X_S, f(X_S), '*g','LineWidth', 2);
        % Точка минимума
        plot(ax, X0, F0, '*r','LineWidth', 4);
        title(ax, sprintf("Точность Eps = %2.0e", e));
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
    arguments
        a    double            % Левая граница отрезка
        b    double            % Правая граница отрезка
        f    function_handle   % Целевая функция
        eps  double            % Точность
    end
    
    x0 = a; x_prev = a; x_next = b;
    f_prev = f(x_prev); f_next = f(x_next); f_m = f(x0);
    d2ff = d2f(f_prev, f_m, f_next, abs(x_next - x_prev));
    
    % f1 = f(a); f2 = f(b); 
    X0 = [x0];
    n = 3; % Число обращений к целевой функции
    
    while true
        dff = df(f_prev, f_next, abs(x_next - x_prev));
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

function [X, F, A_S, B_S, N] = GoldenRatio(a, b, f, eps)
    arguments
        a   double           % Левая граница отрезка
        b   double           % Левая граница отрезка
        f   function_handle  % Целевая функция
        eps double           % Точность
    end
    
    tau = 0.61803;

    A_S = [];
    B_S = [];

    l = b - a;

    x1 = b - tau * l; 
    x2 = a + tau * l; 

    f1 = f(x1);
    f2 = f(x2);
    N = 2;

    while true
        A_S = [A_S, a];
        B_S = [B_S, b];

        if (f1 >= f2)
            a = x1;
            l = b - a;

            x1 = x2; 
            f1 = f2;

            x2 = a + tau * l;
            f2 = f(x2);
        else
            b = x2;
            l = b - a;

            x2 = x1;
            f2 = f1;

            x1 = b - tau * l;
            f1 = f(x1);
        end
        
        N = N + 1;

        if (l < 2 * eps)
            break;
        end
    end

    X = (a + b) / 2.0;
    F = f(X);
    N = N + 1;
end

% Метод парабол
function [X, F, A_S, B_S, N] =  ParabolaMethod(x1, x2, x3, f, eps)
    arguments
        x1    double           % Первая начальная точка (левая граница отрезка)
        x2    double           % Вторая начальная точка (внутренняя точка отрезка)
        x3   double            % Третья начальная точка (правая граница отрезка)
        f    function_handle   % Целевая функция
        eps  double            % Точность
    end
    
    A_S = [];
    B_S = [];
    
    f1 = f(x1); f2 = f(x2); f3 = f(x3);
    x_curr = x2;
    f_s = f2;
    N = 3; % Число обращений к целевой функции
    
    % Повторять, пока отрезок или разность между старым
    % и новым значением f(x*) не станет меньше точности
    while true
        A_S = [A_S, x1];
        B_S = [B_S, x3];
        
        r12 = x1^2 - x2^2;
        r23 = x2^2 - x3^2;
        r31 = x3^2 - x1^2;
        
        s12 = x1 - x2;
        s23 = x2 - x3;
        s31 = x3 - x1;
        
        % Вычисление точки минимума аппроксимирующей параболы
        x_prev = x_curr;
        x_curr = 0.5 * (f1*r23 + f2*r31 + f3*r12) / (f1*s23 + f2*s31 + f3*s12);
        f_s = f(x_curr); N = N + 1;
        
        % Проверка условия окончания поиска
        if (x3 - x1) <= eps || abs(x_prev - x_curr) <= eps
            break;
        end
        
        % Выбор тройки точек для следующей итерации
        if x_curr >= x2 && x_curr <= x3
            if f_s <= f2
                x1 = x2; x2 = x_curr;
                f1 = f2; f2 = f_s;
            else
                x3 = x_curr;
                f3 = f_s;
            end
        elseif x_curr >= x1 && x_curr <= x2
            if f_s <= f2
                x3 = x2; x2 = x_curr;
                f3 = f2; f2 = f_s;
            else
                x1 = x_curr;
                f1 = f_s;
            end
        end
    end
    
    X = x_curr;
    F = f_s;
end

% Метод золотого сечения для нахождения начальных точек метода парабол
function [X1, X2, X3, A_S, B_S, N] =  GoldenRatioForParabola(a, b, f)
    arguments
        a   double           % Левая граница отрезка
        b   double           % Левая граница отрезка
        f   function_handle  % Целевая функция
    end
    
    tau = (5^0.5 - 1) / 2;

    A_S = [];
    B_S = [];

    l = b - a;

    % Вычисление пробных точек
    x1 = b - tau * l; 
    x2 = a + tau * l; 

    f1 = f(x1);
    f2 = f(x2);
    
    N = 2;

    while true
        A_S = [A_S, a];
        B_S = [B_S, b];
        
        % Выбор новых точек для следующей итерации
        if (f1 >= f2)
            a = x1;
            l = b - a;

            x1 = x2; 
            f1 = f2;

            x2 = a + tau * l;
            f2 = f(x2);
            N = N + 1;
        else
            b = x2;
            l = b - a;

            x2 = x1;
            f2 = f1;

            x1 = b - tau * l;
            f1 = f(x1);
            N = N + 1;
        end
        
        % Проверка условия выбора начальных точек для метода парабол
        if (x1 < x2) && (x2 < b) && (f2 <= f1) && (f2 <= f(b))
            X1 = x1; 
            X2 = x2; 
            X3 = b;
            break;
        elseif (a < x1) && (x1 < x2) && (f1 <= f(a)) && (f1 <= f2)
            X1 = a; 
            X2 = x1; 
            X3 = x2;
            break;
        end
     end
end