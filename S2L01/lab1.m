function lab1()
    f = @(x) (tan((x.^4.0 + 2.0.*x.^2.0 - 2.0.*x + 2.0^(1.0/2.0) + 1.0) ./ 8.0) + ...
          sin((4.0.*x.^3.0 - 7.0.*x - 9.0) ./ (20.0.*x + 28.0)));
      
    a = 0; b = 1;
    X = a:0.01:b;
    Y = f(X);
    eps = [0.01, 0.0001, 0.000001];
    
    fprintf("\n���������� ���������� ����� ��������\n");
    fprintf("��� ��������� �������� ��������:\n\n");
    fprintf("Eps - ��������\n");
    fprintf("N - ����� ��������� � ������� �������\n");
    fprintf("x* - ��������� ����� �������� �������\n");
    fprintf("f(x*) - ��������� ������� �������\n\n");
    fprintf(" # |    Eps   | N  |    x*   |   f(x*)\n");
    fprintf("---|----------|----|---------|--------\n");
    
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    title('����� ������������ ������');
    tiledlayout(2, 2);
    
    for i = 1:length(eps)
        % ���������� ����� �������� � �������� �������
        [X0, F0, X_S, F_S, N] = BitwiseSearch(0, 1, f, eps(i));
        
        % ����� ������ ������� ����������� ����������
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, eps(i), N, X0, F0);
        
        % ����� ������� ��� ������ ��������
        ax = nexttile;
        % ������� �������
        plot(ax, X, Y, '-b','LineWidth',1.5);
        hold on;
        % ������������������ �����������
        plot(ax, X_S, F_S, '*g','LineWidth', 2);
        % ����� ��������
        plot(ax, X0, F0, '*r','LineWidth', 4);
        title(ax, sprintf("�������� Eps = %2.0e", eps(i)));
        legend('������� �������', '������������������ �����������', ...
            '����� ��������');
    end
    
end

% ����� ������������ ������
function [X, F, X_S, F_S, N] = BitwiseSearch(a, b, f, eps)
    arguments
        a   double           % ����� ������� �������
        b   double           % ����� ������� �������
        f   function_handle  % ������� �������
        eps double           % ��������
    end
    
    delta = (b - a) / 4.0;  % ��� ������

    x0 = a;      % ��������� ����������� ����� ��������
    f0 = f(x0);  % ��������� ����������� �������� �������
    N = 1;
    
    X_S = [];    % ������ ���� ����������� ����� ��������
    F_S = [];    % ������ ���� ����������� �������� �������

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

  
