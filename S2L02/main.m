method = GoldenRatio();

f = @(x) (sin((x.^4 + x.^3 - 3 * x + 3 - 30.^(1.0/3.0)) / 2.0) + ...
          tanh((4 * 3.^0.5 * x.^3 - 2 * x - 6 * 2.^0.5 + 1) /  ...
               (-2 * 3.^0.5 * x.^3 + x + 3 * 2.^0.5)) + 1.2 );

X = 0:0.01:1;
Y = f(X);
[X0, F0, A_S, B_S] = method.Solve(0, 1, f);

X0
F0
N = length(A_S)

Ranges = [A_S; B_S]

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
title('Метод золотого сечения');
hold on;

plot(X, Y, '-b','LineWidth',1.5);
plot(X0, F0, '*r','LineWidth', 5);
for i = 1:length(A_S)
    plot(A_S(i), f(A_S(i)), '*','LineWidth', 2);
    plot(B_S(i), f(B_S(i)), '*','LineWidth', 2);
end
legend('Целевая функция', 'Точка минимума');
