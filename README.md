# SparseGaussElimination
This Repo is an efficient Fortran90 implementation for directly solving a set of $n$ linear equations $A\mathbf{x}=\mathbf{b}$, where $A\in\mathbb{R}^{n\times n}$ is a $n\times n$ sparse square matrix, $\mathbf{x}\in\mathbb{R}^n$ the solution vector with $n$ elements and $b\in\mathbb{R}^n$ the right hand side of linear equation system. The method is using Gaussian elimination to decompose the matrix $A$ into lower and upper triangular matrices $L$ and $U$, where $L\cdot U=A$. 

Solving very large sets of linear equations is of order $\mathbb{O}(n^3)$. This requires clever usage of the matrix features, e.g. sparsity of the matrix. In several applications (mathematical models of real world problems) large sparse matrices appear, where the number of non-zeros in the matrix $A$ is on the order of one percent or lower. 

This implementation exploits the matrix non-zero structure, reorders rows and columns to avoid fill-in during factorization with a low-cost greedy heuristic introduced by [Markowitz 1957](http://www.jstor.org/stable/2627454), where the most resource consuming part of the method (i.e. obtaining the non-zero pattern of the factors $L$ and $U$) is done only once during initialization. 

**The following steps have to be concidered:**
0. Symbolic Decomposition
  - computes the non-zero pattern of the triangular matrices $L$ and $U$
  - minimum degree heuristic by Markowitz 1957 is applied to avoid massive fill in
1. Numerical Factorization
  - numbers are calculated 
3. Foward-Backward-Solve
  - solves lower and upper triangular subsystems to obtain the solution vector $x$


