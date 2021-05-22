function lab3()
    f = @(x) (tan((x.^4 + 2.*x.^2 - 2.*x + 2^0.5 + 1) ./ 8) + ...
          sin((4.*x.^3 - 7.*x - 9) ./ (20.*x + 28)));
      
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
        [A, B, C, A_S, B_S] = GoldenRatio(0, 1, f, 0.01);
        [X0, F0, A_S, B_S] = ParabolaMethod(A, B, C, f, eps(i));
        
        % ����� ������ ������� ����������� ����������
        fprintf("%2i | %5f | %2i | %5.5f | %5.5f\n", ...
            i, eps(i), length(A_S) + 1, X0, F0);
        
        % ����� ������� ��� ������ ��������
        ax = nexttile;
        % ������� �������
        plot(ax, X, Y, '-b','LineWidth',1.5);
        hold on;
        % ������������������ �����������
        plot(ax, [A_S, B_S], [f(A_S), f(B_S)], '*g','LineWidth', 2);
        % ����� ��������
        plot(ax, X0, F0, '*r','LineWidth', 4);
        title(ax, sprintf("�������� Eps = %2.0e", eps(i)));
        legend('������� �������', '������������������ �����������', ...
            '����� ��������');
    end
    
end

function [X, F, A_S, B_S] =  ParabolaMethod(a, b, c, f, eps)
    arguments
        a    double
        b    double
        c   double           % ����� ������� �������
        f    function_handle   % ������� �������
        eps  double            % ��������
    end
    
    A_S = [];
    B_S = [];
    
    x1 = a; x2 = c; x3 = b;
    
    f1 = f(x1); f2 = f(x2); f3 = f(x3);
    n = 3; % ����� ��������� � ������� �������
    
    % ���������, ���� ������� ��� �������� ����� ������
    % � ����� ��������� f(x*) �� ������ ������ ��������
    while (x3 - x1) > eps
        A_S = [A_S, x1];
        B_S = [B_S, x3];
        
        r12 = x1^2 - x2^2;
        r23 = x2^2 - x3^2;
        r31 = x3^2 - x1^2;
        
        s12 = x1 - x2;
        s23 = x2 - x3;
        s31 = x3 - x1;
        
        % ���������� ����� �������� ���������������� ��������
        x_s = 0.5 * (f1*r23 + f2*r31 + f3*r12) / (f1*s23 + f2*s31 + f3*s12);
        f_s = f(x_s); n = n + 1;
        
        % ����� ������ ����� ��� ��������� ��������
        if x_s >= x2 && x_s <= x3
            if f_s <= f2
                x1 = x2; x2 = x_s;
                f1 = f2; f2 = f_s;
            else
                x3 = x_s;
                f3 = f_s;
            end
        elseif x_s >= x1 && x_s <= x2
            if f_s <= f2
                x3 = x2; x2 = x_s;
                f3 = f2; f2 = f_s;
            else
                x1 = x_s;
                f1 = f_s;
            end
        end
    end
    
    X = x_s;
    F = f_s;
end

function [A, B, C, A_S, B_S] =  GoldenRatio(a, b, f, eps)
    arguments
        a   double           % ����� ������� �������
        b   double           % ����� ������� �������
        f   function_handle  % ������� �������
        eps double           % ��������
    end
    
    tau = (5^0.5 - 1) / 2;

    A_S = [];
    B_S = [];

    l = b - a;

    x1 = b - tau * l; 
    x2 = a + tau * l; 

    f1 = f(x1);
    f2 = f(x2);

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

        if (l < 2 * eps)
            break;
        end
    end
    
    A = a;
    B = b;
    C = (a + b) / 2;
end
