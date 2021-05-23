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
    title('Метод парабол');

    for i = 1:length(eps)
        % Вычисление точки минимума и минимума функции
        [X0, F0, A_S, B_S, N] = ParabolaMethod(a, b, f, eps(i));
        
        % Вывод строки таблицы результатов вычислений
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, eps(i), N, X0, F0);
        
        % Вывод графика для данной точности
        subplot(2, 2, i);
        % Целевая функция
        plot(X, Y, '-b','LineWidth',1.5);
        hold on;
        % Последовательность приближений
        plot([A_S, B_S], [f(A_S), f(B_S)], '*g','LineWidth', 2);
        % Точка минимума
        plot(X0, F0, '*r','LineWidth', 4);
        title(sprintf("Точность Eps = %2.0e", eps(i)));
        legend('Целевая функция', 'Последовательность приближений', ...
            'Точка минимума');
    end
    
end

% Метод парабол
function [X, F, A_S, B_S, N] =  ParabolaMethod(a, b, f, eps)
    [X0, F0, A_S, B_S, N] = GoldenRatio(a, b, f);
    
    [x1, x2, x3] = deal(X0(1), X0(2), X0(3));
    [f1, f2, f3] = deal(F0(1), F0(2), F0(3));
    
    x_curr = x2;
    f_s = f2;
    is_first_iter = true;
    
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
        x_curr = 0.5 * (f1*r23 + f2*r31 + f3*r12) / ...
            (f1*s23 + f2*s31 + f3*s12);
        
        % Проверка условия окончания поиска (на всех итерациях кроме первой)
        if is_first_iter
            is_first_iter = false;
        elseif abs(x_prev - x_curr) <= eps
            break;
        end
        
        f_s = f(x_curr); 
        N = N + 1;
        
        % Выбор тройки точек для следующей итерации
        if x_curr >= x2 && x_curr <= x3
            if f_s <= f2
                x1 = x2; f1 = f2;
                x2 = x_curr; f2 = f_s;
            else
                x3 = x_curr;
                f3 = f_s;
            end
        elseif x_curr >= x1 && x_curr <= x2
            if f_s <= f2
                x3 = x2; f3 = f2; 
                x2 = x_curr; f2 = f_s;
            else
                x1 = x_curr;
                f1 = f_s;
            end
        end
    end
    
    X = x_curr;
    F = f(x_curr);
    N = N + 1;
end

% Метод золотого сечения для выбора начальных точек метода парабол
function [X, F,  A_S, B_S, N] =  GoldenRatio(a, b, f)
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
        
        % Проверка условия выбора начальных точек для метода парабол
        if (x1 < x2) && (x2 < b) && (f2 <= f1) && (f2 <= f(b))
            X = [x1, x2, b];
            F = [f1, f2, f(b)];
            N = N + 1;
            break;
        elseif (a < x1) && (x1 < x2) && (f1 <= f(a)) && (f1 <= f2)
            X = [a, x1, x2];
            F = [f(a), f1, f2];
            N = N + 1;
            break;
        end
        
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
     end
end
