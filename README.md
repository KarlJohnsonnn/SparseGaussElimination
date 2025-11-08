# SparseGaussElimination
This Repo is an efficient Fortran90 implementation for directly solving a set of $n$ linear equations $A\mathbf{x}=\mathbf{b}$, where $A\in\mathbb{R}^{n\times n}$ is a $n\times n$ sparse square matrix, $\mathbf{x}\in\mathbb{R}^n$ the solution vector with $n$ elements and $b\in\mathbb{R}^n$ the right hand side of linear equation system. The method is using Gaussian elimination to decompose the matrix $A$ into lower and upper triangular matrices $L$ and $U$, where $L\cdot U=A$. 

Solving very large sets of linear equations is of order $\mathbb{O}(n^3)$. This requires clever usage of the matrix features, e.g. sparsity of the matrix. In several applications (mathematical models of real world problems) large sparse matrices appear, where the number of non-zeros in the matrix $A$ is on the order of one percent or lower. 

This implementation exploits the matrix non-zero structure, reorders rows and columns to avoid fill-in during factorization with a low-cost greedy heuristic introduced by [Markowitz 1957](http://www.jstor.org/stable/2627454), where the most resource consuming part of the method (i.e. obtaining the non-zero pattern of the factors $L$ and $U$) is done only once during initialization. 

**The following steps have to be considered:**

1. Symbolic Decomposition
    - computes the non-zero pattern of the triangular matrices $L$ and $U$
    - minimum degree heuristic by Markowitz 1957 is applied to avoid massive fill-in
2. Numerical Factorization
    - numbers are calculated 
3. Forward-Backward-Solve
    - solves lower and upper triangular subsystems to obtain the solution vector $x$

## Compilation

```bash
bash Compile.sh
```

Or manually:
```bash
gfortran -o SparseGaussElimination.exe Kind_Mod.f90 mo_unirnk.f90 Sparse_Mod.f90 Main_SparseGaussElimination.f90
```

## Usage

Run the executable and provide a sparse matrix file when prompted:
```bash
./SparseGaussElimination.exe [optional_matrix_file]
```

Example matrix files are located in the `MATRICES/` directory.

## API Reference

The `sparse_mod` module provides comprehensive functionality for sparse matrix operations. Below is a complete reference of all available functions and subroutines.

### Data Types

| Type | Description |
|------|-------------|
| `csr_matrix_t` | Compressed Sparse Row (CSR) matrix format (main format) |
| `sprowcold_t` | Dynamic sparse row-column format for symbolic factorization |
| `sparse_row_ind_col_ind_t` | Standard row index, column index format |

### Matrix Creation and Management

| Function/Subroutine | Purpose | Signature |
|---------------------|---------|-----------|
| `new_csr` | Create and initialize a new CSR matrix | `function(m, n, [nnz], [ri], [ci], [val])` |
| `sparse_identity` | Create a sparse identity matrix | `function(dim)` |
| `copy_csr` | Create a deep copy of a CSR matrix | `function(orig)` |
| `free_matrix_csr` | Deallocate CSR matrix arrays | `subroutine(a)` |
| `free_sprowcold` | Deallocate dynamic sparse row-column structure | `subroutine(a)` |
| `free_sparse_row_ind_col_ind` | Deallocate sparse row-column index structure | `subroutine(a)` |

### Format Conversions

| Function/Subroutine | Purpose | Signature |
|---------------------|---------|-----------|
| `csr_to_full` | Convert CSR matrix to full dense format | `function(csr)` |
| `full_to_csr` | Convert full dense matrix to CSR format | `function(full)` |
| `csr_to_sprowcold` | Convert CSR to dynamic row-column format | `function(csr)` |
| `rowcold_to_csr` | Convert dynamic row-column format to CSR | `function(sp_row)` |

### Matrix Operations

| Function/Subroutine | Purpose | Signature |
|---------------------|---------|-----------|
| `transpose_sparse` | Compute transpose of a sparse matrix | `subroutine(mat_at, mat_a)` |
| `symbolic_mult` | Symbolic matrix multiplication C = A × B | `subroutine(a, b, c)` |
| `sparse_mult` | Numeric matrix multiplication C = A × B | `subroutine(a, b, c)` |
| `symbolic_add` | Symbolic matrix addition C = A + B | `subroutine(mat_c, mat_a, mat_b)` |
| `sparse_add` | Numeric matrix addition/subtraction C = A ± B | `subroutine(mat_c, mat_a, mat_b, [sub])` |
| `dax_sparse` | Sparse matrix-vector product y = A × x | `function(a, x)` (overloads `*` operator) |
| `daxpy_sparse` | Sparse matrix-vector product with addition y = A × x + y | `function(a, x, y)` |

### LU Factorization and Solving

| Function/Subroutine | Purpose | Signature |
|---------------------|---------|-----------|
| `symblu_sprowcold` | Symbolic LU factorization with given pivot ordering | `subroutine(a, permu)` |
| `symblu_sprowcold_m` | Symbolic LU factorization with Markowitz minimum-degree ordering | `subroutine(a)` |
| `sparse_lu` | Numeric LU factorization (in-place) | `subroutine(a)` |
| `solve_sparse` | Solve LU × x = b using forward/backward substitution | `subroutine(lu, rhs)` |
| `get_lu_permutation` | Get permutation vector mapping original to LU structure | `subroutine(permutation, lu, a)` |

### Input/Output Operations

| Function/Subroutine | Purpose | Signature |
|---------------------|---------|-----------|
| `read_sparse_matrix` | Read sparse matrix from file | `function(filename)` |
| `write_sparse_matrix` | Write sparse matrix to file | `subroutine(a, [filename], [nr], [ns], [extended], [classic], [teq], [linalg])` |
| `print_sparse_matrix` | Print sparse matrix to stdout | `subroutine(a)` |
| `input_pivot_order` | Read pivot ordering from file | `function(filename, dim)` |

### Utility Functions (Internal)

| Function/Subroutine | Purpose | Signature |
|---------------------|---------|-----------|
| `sort_vec` | Simple bubble sort for integer vectors | `subroutine(vec)` |
| `swap_int` | Swap two integer values | `subroutine(i, j)` |
| `insert_sprowcold` | Insert element into dynamic sparse structure | `subroutine(a, i_a, j_a, ins)` |
| `gcmat_sprowcold` | Garbage collection for dynamic sparse structure | `subroutine(a)` |

### Usage Example

```fortran
use kind_mod
use sparse_mod

! Create a matrix
type(csr_matrix_t) :: A, LU
real(dp), allocatable :: rhs(:)

! Read matrix from file
A = read_sparse_matrix('MATRICES/Miter_SmallStratoKPP_cl.SparseMat')

! Convert to dynamic format for factorization
type(sprowcold_t) :: temp_lu
temp_lu = csr_to_sprowcold(A)

! Perform symbolic LU factorization with Markowitz ordering
call symblu_sprowcold_m(temp_lu)

! Convert back to CSR format
LU = rowcold_to_csr(temp_lu)

! Map values from original matrix to LU structure
integer, allocatable :: lu_perm(:)
call get_lu_permutation(lu_perm, LU, A)

! Numeric factorization
call sparse_lu(LU)

! Solve system
allocate(rhs(LU%m))
rhs = 1.0_dp  ! example right-hand side
call solve_sparse(LU, rhs)
! rhs now contains the solution

! Clean up
call free_matrix_csr(A)
call free_matrix_csr(LU)
deallocate(rhs)
```

### Matrix File Format

The sparse matrix file format is simple and human-readable:

```
m n [comment]
nnz [comment]
[blank line]
i j value
i j value
...
```

Where:
- `m` = number of rows
- `n` = number of columns
- `nnz` = number of nonzero entries
- `i j value` = row index, column index, and value for each nonzero entry

## References

- Markowitz, H. M. (1957). "The Elimination Form of the Inverse and its Application to Linear Programming". *Management Science*, 3(3), 255-269. [DOI](http://www.jstor.org/stable/2627454)


