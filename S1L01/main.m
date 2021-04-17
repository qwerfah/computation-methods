A1 = [7, 4, 3, 8, 2;
      4, 5, 1, 6, 3;
      8, 4, 5, 7, 2;
      1, 2, 4, 7, 2;
      3, 9, 9, 2, 5];

A = [7, 7, 4, 6, 5; 
     3, 8, 1, 8, 8; 
     5, 5, 7, 4, 1; 
     7, 6, 8, 6, 3; 
     4, 9, 2, 4, 3];
 
method = HungarianMethod(A);
[srow, scol] = method.BuildZeroSystem(true);
[X, f] = method.CalcOptSolution(srow, scol);

X
f
 




