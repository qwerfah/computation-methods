classdef MatrixOperations
    methods(Static)
        % Вычесть минимумы в каждой строке матрицы
        function Result = SubtractMinInRow(matrix)
            arguments
                matrix (:,:) double
            end
            
            tr = transpose(matrix);
            
            Result = transpose(tr - min(tr));
        end
        % Вычесть минимумы в каждом столбце матрицы
        function Result = SubtractMinInColumn(matrix)
            arguments
                matrix (:,:) double
            end
            
            Result = matrix - min(matrix);
        end
        
        % Построить систему независимых нулей
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
        
        % Улучшить систему независимых нулей
        function [Matrix, sRow, sCol] = UpgradeZeroSystem(matrix, row, col)
            arguments
                matrix (:,:) double
                row, col (1,:) double
            end
            
            % Столбцы и строки с 0*
            mrow = []; mcol = col;
            % Индексы элементов, помеченных '
            srow = []; scol = [];
            
            fprintf("\nОтмечаем столбцы, содержащие 0*:\n");
            MatrixOperations.PrintMatrix(...
                matrix, row, col, srow, scol, mrow, mcol, [], []);
            
            while true
                % Находим среди неотмеченных элементов 0
                [found, x, y] = MatrixOperations.FindNonMarkedZero(...
                    matrix, mrow, mcol);
                % Если не найден, вычитаем из невыделенных столбцов и
                % прибавляем к выделенным строкам минимум среди
                % невыделенных элементов
                while ~found
                    matrix = MatrixOperations.SubtractNonMarkedMin(...
                        matrix, mrow, mcol);
                    
                    fprintf("\nНет нулей, вычитаем минимум:\n");
                    MatrixOperations.PrintMatrix(...
                        matrix, row, col, srow, scol, mrow, mcol, [], []);
                    
                    [found, x, y] = MatrixOperations.FindNonMarkedZero(...
                        matrix, mrow, mcol);
                end
                % Добавляем новый 0' в список помеченных '
                srow = [srow, x]; scol = [scol, y];
                
                fprintf("\nОтмечаем невыделенный 0 символом ':\n");
                MatrixOperations.PrintMatrix(...
                    matrix, row, col, srow, scol, mrow, mcol, [], []);
                
                % Если в одной строке с 0' есть 0*
                if any(row == x)
                    ind = mcol == col(row == x);
                    % Снять выделение со столбца с 0*
                    mcol(ind) = [];
                    % Пометить строку с 0'
                    mrow = [mrow, x];
                    
                    fprintf("\nВ строке с 0' есть 0*, снимаем " + ...
                        "выделение со столбца, выделяем строку:\n");
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
        
        % Найти непомеченный 0 среди невыделенных элементов матрицы
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
        
        % Вычесть из невыделенных столбцов и прибавить к выделенным 
        % строкам минимум среди невыделенных элементов
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
        
        % Построить непродолжаемую L-цепочку и обновить систему нулей
        function [sRow, sCol] = BuildLChain(matrix, row, col, srow, scol)
            arguments
                matrix (:,:) double
                row, col (1,:) double
                srow, scol (1,:) double
            end
            
            % Начиная с последнего отмеченного 0'
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
            
            fprintf("\nВ строке с 0' нет 0*, строим L-цепочку:\n");
            MatrixOperations.PrintMatrix(...
                matrix, row, col, srow, scol, [], [], lrow, lcol);
            
            for i = 1:length(lrow)
                % Если 0*, то снимаем *
                if mod(i, 2) == 0 
                    rowind = row == lrow(i);
                    colind = col == lcol(i);
                    row(colind & rowind) = [];
                    col(colind & rowind) = [];
                % Если 0', то помечаем *
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