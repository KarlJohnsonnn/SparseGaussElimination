# SparseGaussElimination
This Repo is a highly efficient Fortran90 implementation to directly solve a set of $N$ linear equations $A\mathbr{x}=\mathbr{b}$, where $A\in\mathcal{R}^{n\times n}$ is a $n\times n$ sparse square matrix, $\mathbr{x}\in\mathcal{R}^n$ the solution vector with $n$ elements and $b\in\mathcal{R}^n$ the right hand side of linear equation system. The method is using Gaussian elimination to decompose the matrix $A$ into lower and upper triangular matrices $L$ and $U$, where $L\cdot U=A$. 

Solving very large sets of linear equations of order $\mathcal{O}(n^3)$ can be computationally very expensive and requires fast implementations and clever usage of matrix features, e.g. sparsity of the matrix. In several applications (mathematical models) large sparse matrices appear, where the number of non-zeros in the matrix $A$ is on the order of one percent (or even lower). 

This implementation exploits the matrix non-zero structure, reorders rows and columns to avoid fill-in during factorization using a low-cost greedy heuristic introduced by [Markowitz 1957](http://www.jstor.org/stable/2627454). The expensive part of the method (getting the non-zero pattern of the decomposition) is done only once during initialization. 

**The following steps have to be concidered:**
0. Symbolic Decomposition
  - computes the non-zero pattern of the triangular matrices $L$ and $U$
  - minimum degree heuristic by Markowitz 1957 is applied to avoid massive fill in
1. Numerical Factorization
  - numbers are calculated 
3. Foward-Backward-Solve
  - solves lower and upper triangular subsystems to obtain the solution vector $x$


