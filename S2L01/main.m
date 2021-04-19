method = BitwiseSearch();

f = @(x) (sin((x.^4.0 + x.^3.0 - 3.0 * x + 3.0 - 30.^(1.0/3.0)) / 2.0) + ...
          tanh((4.0 * 3.0.^0.5 * x.^3.0 - 2.0 * x - 6.0 * 2.0.^0.5 + 1.0) /  ...
               (-2.0 * 3.0.^0.5 * x.^3.0 + x + 3.0 * 2.0.^0.5)) + 1.2 );

X = 0:0.01:1;
Y = f(X);
[X0, F0, X_S, F_S] = method.Solve(0, 1, f);

X0
F0
N = length(X_S)

i = 5;
F_S - f(X_S)

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
title('Метод поразрядного поиска');
hold on;

plot(X, Y, '-b','LineWidth',1.5);
plot(X0, F0, 'or','LineWidth', 5);
plot(X_S, F_S, 'og','LineWidth', 2);
legend('Целевая функция', 'Точка минимума', 'Последовательность приближений');
