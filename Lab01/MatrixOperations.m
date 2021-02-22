classdef MatrixOperations
    methods(Static)
        % ������� �������� � ������ ������ �������
        function Result = SubtractMinInRow(matrix)
            arguments
                matrix (:,:) double
            end
            
            tr = transpose(matrix);
            
            Result = transpose(tr - min(tr));
        end
        % ������� �������� � ������ ������� �������
        function Result = SubtractMinInColumn(matrix)
            arguments
                matrix (:,:) double
            end
            
            Result = matrix - min(matrix);
        end
        
        % ��������� ������� ����������� �����
        function [sRow, sCol] = BuildZeroSystem(matrix)
            arguments
                matrix (:,:) double
            end
            
            srow = []; scol = [];
            [row, col] = find(~matrix);
            
            for i = 1:length(row)
                if (~any(scol == col(i))) && (~any(srow == row(i)))
                    srow = [srow, row(i)];
                    scol = [scol, col(i)];
                end
            end
            
            sRow = srow; sCol = scol;
        end
        
        % �������� ������� ����������� �����
        function [Matrix, sRow, sCol] = UpgradeZeroSystem(matrix, row, col)
            arguments
                matrix (:,:) double
                row, col (1,:) double
            end
            
            % ������� � ������ � 0*
            mrow = []; mcol = col;
            % ������� ���������, ���������� '
            srow = []; scol = [];
            
            fprintf("\n�������� �������, ���������� 0*:\n");
            MatrixOperations.PrintMatrix(...
                matrix, row, col, srow, scol, mrow, mcol, [], []);
            
            while true
                % ������� ����� ������������ ��������� 0
                [found, x, y] = MatrixOperations.FindNonMarkedZero(...
                    matrix, mrow, mcol);
                % ���� �� ������, �������� �� ������������ �������� �
                % ���������� � ���������� ������� ������� �����
                % ������������ ���������
                while ~found
                    matrix = MatrixOperations.SubtractNonMarkedMin(...
                        matrix, mrow, mcol);
                    
                    fprintf("\n��� �����, �������� �������:\n");
                    MatrixOperations.PrintMatrix(...
                        matrix, row, col, srow, scol, mrow, mcol, [], []);
                    
                    [found, x, y] = MatrixOperations.FindNonMarkedZero(...
                        matrix, mrow, mcol);
                end
                % ��������� ����� 0' � ������ ���������� '
                srow = [srow, x]; scol = [scol, y];
                
                fprintf("\n�������� ������������ 0 �������� ':\n");
                MatrixOperations.PrintMatrix(...
                    matrix, row, col, srow, scol, mrow, mcol, [], []);
                
                % ���� � ����� ������ � 0' ���� 0*
                if any(row == x)
                    ind = mcol == col(row == x);
                    % ����� ��������� �� ������� � 0*
                    mcol(ind) = [];
                    % �������� ������ � 0'
                    mrow = [mrow, x];
                    
                    fprintf("\n� ������ � 0' ���� 0*, ������� " + ...
                        "��������� �� �������, �������� ������:\n");
                    MatrixOperations.PrintMatrix(...
                        matrix, row, col, srow, scol, mrow, mcol, [], []);
                else
                    [sRow, sCol] = MatrixOperations.BuildLChain(...
                        matrix, row, col, srow, scol);
                    Matrix = matrix;
                    break;
                end
            end
        end
        
        % ����� ������������ 0 ����� ������������ ��������� �������
        function [found, x, y] = FindNonMarkedZero(matrix, row, col)
            arguments
                matrix (:,:) double
                row, col (1,:) double
            end
            
            for i = setdiff(1:size(matrix, 1), row)
                for j = setdiff(1:size(matrix, 2), col)
                    if matrix(i, j) == 0
                        x = i; y = j; found = true;
                        return;
                    end
                end
            end
            
            [found, x, y] = deal(false, -1, -1);
        end
        
        % ������� �� ������������ �������� � ��������� � ���������� 
        % ������� ������� ����� ������������ ���������
        function Result = SubtractNonMarkedMin(matrix, row, col)
            arguments
                matrix (:,:) double
                row, col (1,:) double
            end
            
            nrow = setdiff(1:size(matrix, 1), row);
            ncol = setdiff(1:size(matrix, 2), col);
            min = max(matrix, [], 'all');
            
            for i = nrow
                for j = ncol
                    if matrix(i, j) < min
                        min = matrix(i, j);
                    end
                end
            end
            
            for i = ncol
                for j = 1:size(matrix, 1)
                    matrix(j, i) = matrix(j, i) - min;
                end
            end
            
            for i = row
                for j = 1:size(matrix, 2)
                    matrix(i, j) = matrix(i, j) + min;
                end
            end
            
            Result = matrix;
        end
        
        % ��������� �������������� L-������� � �������� ������� �����
        function [sRow, sCol] = BuildLChain(matrix, row, col, srow, scol)
            arguments
                matrix (:,:) double
                row, col (1,:) double
                srow, scol (1,:) double
            end
            
            % ������� � ���������� ����������� 0'
            [x, y] = deal(srow(length(srow)), scol(length(scol)));
            lrow = x; lcol = y;
            
            while true
                x = row(col == y);
                if isempty(x)
                    break;
                end
                lrow = [lrow, x]; lcol = [lcol, y];
                y = scol(srow == x);
                lrow = [lrow, x]; lcol = [lcol, y];
            end
            
            fprintf("\n� ������ � 0' ��� 0*, ������ L-�������:\n");
            MatrixOperations.PrintMatrix(...
                matrix, row, col, srow, scol, [], [], lrow, lcol);
            
            for i = 1:length(lrow)
                % ���� 0*, �� ������� *
                if mod(i, 2) == 0 
                    rowind = row == lrow(i);
                    colind = col == lcol(i);
                    row(colind & rowind) = [];
                    col(colind & rowind) = [];
                % ���� 0', �� �������� *
                else 
                    row = [row, lrow(i)];
                    col = [col, lcol(i)];
                end
            end
            
            [sRow, sCol] = deal(row, col);
        end
        
        function PrintMatrix(...
                matrix, row, col, srow, scol, mrow, mcol, lrow, lcol)
            arguments
                matrix (:,:) double
                row, col (1,:) double
                srow, scol (1,:) double
                mrow, mcol (1,:) double
                lrow, lcol (1,:) double
            end
            
            for i = 1:size(matrix, 1)
                if any(mcol == i)
                    fprintf('    + ');
                else
                    fprintf('      ');
                end
            end
            fprintf('\n');
            
            for i = 1:size(matrix, 1)
                for j = 1:size(matrix, 2)
                    if any(lrow == i & lcol == j) && ...
                       any(row == i & col == j)
                        fprintf('%3d*L ', matrix(i, j));
                    elseif any(lrow == i & lcol == j) && ...
                           any(srow == i & scol == j)
                        fprintf('%3d''L ', matrix(i, j));
                    elseif any(row == i & col == j)
                        fprintf('%4d* ', matrix(i, j));
                    elseif any(srow == i & scol == j)
                        fprintf('%4d'' ', matrix(i, j));
                    else
                        fprintf('%5d ', matrix(i, j));
                    end 
                end
                if any(mrow == i)
                    fprintf(' +\n');
                else
                    fprintf('\n');
                end
            end
        end
    end
end