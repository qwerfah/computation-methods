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

% ����� ������� ��� ���������� ������� fminbnd
function FminbndTable(f, a, b, eps)
    fprintf("\n���������� ���������� ����� ��������\n");
    fprintf("��� ��������� �������� �������� �������� fminbnd:\n\n");
    fprintf(" # |    Eps   | N  |    x*   |   f(x*)\n");
    fprintf("---|----------|----|---------|--------\n");
    for i = 1:length(eps)
        options = optimset('TolX', eps(i));
        [x, fval, ~, output] = fminbnd(f, a, b, options);
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, eps(i), output.funcCount, x, fval);
    end
end

% ����� ������� � �������� ��� ����������� ������ ������ �������
function NeutonMethodTable(f, a, b, X, Y, eps)
    fprintf("\n���������� ���������� ����� ��������\n");
    fprintf("��� ��������� �������� �������� ������� �������:\n\n");
    fprintf("Eps - ��������\n");
    fprintf("N - ����� ��������� � ������� �������\n");
    fprintf("x* - ��������� ����� �������� �������\n");
    fprintf("f(x*) - ��������� ������� �������\n\n");
    fprintf(" # |    Eps   | N  |    x*   |   f(x*)\n");
    fprintf("---|----------|----|---------|--------\n");
    
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    title('����� �������');
    
    for i = 1:length(eps)
        % ���������� ����� �������� � �������� �������
        [X0, F0, X_S, N] = NeutonMethod(a, b, f, eps(i));
        
        % ����� ������ ������� ����������� ����������
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, eps(i), N, X0, F0);
        
        % ����� ������� ��� ������ ��������
        subplot(2, 2, i);
        % ������� �������
        plot(X, Y, '-b','LineWidth',1.5);
        hold on;
        % ������������������ �����������
        plot(X_S, f(X_S), '*g','LineWidth', 2);
        % ����� ��������
        plot(X0, F0, '*r','LineWidth', 4);
        title(sprintf("�������� Eps = %2.0e", eps(i)));
        legend('������� �������', '������������������ �����������', ...
            '����� ��������');
    end
end

% ���������� ������������� ������ �����������
function df = df(f_prev, f_next, delta)
    df = (f_next - f_prev) / 2 / delta;
end

% ���������� ������������� ������ �����������
function d2f = d2f(f_prev, f_m, f_next, delta)
    d2f = (f_next - 2*f_m + f_prev) / delta / delta;
end

% ����� �������
function [X, F, X0, N] =  NeutonMethod(a, b, f, eps)
    x0_prev = (a + b) / 1.3;
    x_prev = x0_prev - eps;
    x_next = x0_prev + eps;
    
    f_prev = f(x_prev);
    f_next = f(x_next);
    
    d2ff = d2f(f_prev, f(x0_prev), f_next, x_next - x_prev);

    X0 = [x0_prev];
    N = 3; % ����� ��������� � ������� �������
    
    while true
        % ��������� ����������� ����� ��������
        dff = df(f_prev, f_next, x_next - x_prev);
        x0 = x0_prev - dff / d2ff;
        
        % �������� ������� ���������� ������
        if abs(x0_prev - x0) <= eps
            break;
        end
        
        % ���������� ����� ��� ���������� ������������� �����������
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