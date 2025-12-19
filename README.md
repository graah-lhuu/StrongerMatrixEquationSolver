# StrongerMatrixEquationSolver
A stronger matrix solver which can compute x hat if x cannot be solved.  
When there is no solution, we see whether A is full column rank.  
If yes, we can use a simple formula A^T * A * x_hat = A^T*b.  
