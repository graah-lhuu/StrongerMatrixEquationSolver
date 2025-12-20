# StrongerMatrixEquationSolver
This is a stronger matrix solver, which can compute x_hat if x cannot be solved and A is full column rank.  
When there is no solution, we see whether A is full column rank.  
If yes, we can use a simple formula A^T * A * x_hat = A^T * b.  
However, this x_hat is not the previous x.  
The equation can be solved exactly when vector b is in the column space of A.  
When b is not in the col(A), we need to do a projection of b.  
More precisely, we change the b into p, which is a projection of b onto column space.  
(By the way, apparently, the error vector e = b - p, and vector e must be in the null space of A transpose.)  
Actually, A * x_hat = p definitely can be solved since p is in the column space of A, and x_hat is the solution to A * x_hat = p.  
In this StrongMatrixEquationSolver, I removed the examples and focused on the solutions.
