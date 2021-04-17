classdef HungarianMethod
    properties
        Matrix (:,:) double
    end
    
    methods
        function obj = HungarianMethod(matrix)
            obj.Matrix = matrix;
        end
        
        function [srow, scol] = BuildZeroSystem(obj, isMax)
            if isMax
                matrix = max(obj.Matrix, [], 'all') - obj.Matrix;
            else
                matrix = obj.Matrix;
            end
            
            fprintf("Исходная матрица:\n");
            MatrixOperations.PrintMatrix(...
                matrix, [], [], [], [], [], [], [], []);
            
            matrix = MatrixOperations.SubtractMinInColumn(matrix);
            fprintf("\nВычитание минимумов в столбцах:\n");
            MatrixOperations.PrintMatrix(...
                matrix, [], [], [], [], [], [], [], []);
            
            matrix = MatrixOperations.SubtractMinInRow(matrix);
            fprintf("\nВычитание минимумов в строках:\n");
            MatrixOperations.PrintMatrix(...
                matrix, [], [], [], [], [], [], [], []);
            
            [srow, scol] = MatrixOperations.BuildZeroSystem(matrix);
            fprintf("\nСистема независымых нулей №1:\n");
            MatrixOperations.PrintMatrix(...
                matrix, srow, scol, [], [], [], [], [], []);

            i = 2;
            while length(srow) ~= size(matrix, 1)
                [matrix, srow, scol] = ...
                    MatrixOperations.UpgradeZeroSystem(matrix, srow, scol);
                fprintf("\nСистема независымых нулей №%d:\n", i);
                MatrixOperations.PrintMatrix(...
                    matrix, srow, scol, [], [], [], [], [1,2,3], [3,2,1]);
                i = i + 1;
            end
        end
        
        function [X, f] = CalcOptSolution(obj, row, col)
            arguments
                obj
                row, col (1,:) double
            end
            
            X = zeros(size(obj.Matrix, 1));
            for i = 1:length(row)
                X(row(i), col(i)) = 1;
            end
            
            f = sum(obj.Matrix .* X, 'all');
        end
        
        function Result = EvalDebug(obj)
            Result = 0;
        end
    end
end