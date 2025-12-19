# StrongerMatrixEquationSolver
A stronger matrix solver which can compute x_hat if x cannot be solved.  
When there is no solution, we see whether A is full column rank.  
If yes, we can use a simple formula A^T * A * x_hat = A^T * b.  
However, this x_hat is not the previous x.  
More precisely, we change the b into p, which is a projection of b onto column space.  
Actually, x_hat is the solution to A * x_hat = p.  
In this StrongMatrixEquationSolver, I removed the examples and focused on the solutions.
