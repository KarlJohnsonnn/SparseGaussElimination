!===============================================================================
! Program: Main_SparseGaussElimination
!
! Purpose: Test program for the solution of a sparse linear equation system
!
! Description:
!   Solves Ax = b where A is a sparse matrix using LU decomposition with
!   Markowitz ordering for optimal fill-in reduction during factorization.
!
! Author: Willi Schimmel
! Institute: Leibniz Institute for Tropospheric Research (TROPOS)
! Date: 7 June 2018
!===============================================================================
program main_sparse_gauss_elimination

  use kind_mod
  use sparse_mod

  implicit none
  
  character(len=80) :: filename0 = ''
  character(len=80) :: input = ''
  character(len=400) :: io_msg
  integer :: io_stat
  integer :: i
  
  ! Matrices for symbolic phase
  type(csr_matrix_t) :: miter, lu_miter
  type(sprowcold_t) :: temp_lu_dec
  real(dp), allocatable :: rhs(:)
  integer, allocatable :: lu_perm(:)
  integer, allocatable :: piv_order(:)
  
  ! Timer variables
  real(dp) :: start_timer, time_symbolic, time_numeric, time_solve
  
  character(len=*), parameter :: fmt_output = '(10X,A)'

  ! Read matrix filename from command line or prompt user
  call getarg(1, filename0)
  if (filename0 == '') then
    write(*,fmt_output,advance='no') 'Input matrix filename: '
    read(*,*,iostat=io_stat,iomsg=io_msg) filename0
  end if
  filename0 = trim(adjustl(filename0))
  
  miter = read_sparse_matrix(filename0)

  call cpu_time(start_timer)
  
  ! Convert to sparse row-column format for factorization
  temp_lu_dec = csr_to_sprowcold(miter)

  ! Symbolic factorization - choose ordering method
  write(*,*)
  write(*,fmt_output) 'Ordering options, type:'
  write(*,*)
  write(*,fmt_output) '      marko     ==> Markowitz ordering heuristic (minimum-degree)'
  write(*,fmt_output) '      given     ==> pivot order is provided by an integer array'
  write(*,fmt_output) '      inorder   ==> pivot order is not altered'
  write(*,*)
  write(*,fmt_output,advance='no') 'Choose ordering :: '
  read(*,*,iostat=io_stat,iomsg=io_msg) input

  if (trim(input) == 'marko') then
    call symblu_sprowcold_m(temp_lu_dec)
  
  else if (trim(input) == 'given') then
    write(*,fmt_output,advance='no') 'Select file :: '
    read(*,*,iostat=io_stat,iomsg=io_msg) input
    piv_order = input_pivot_order(input, temp_lu_dec%m)
    call symblu_sprowcold(temp_lu_dec, piv_order)
  
  else if (trim(input) == 'inorder') then
    piv_order = [(i, i = 1, temp_lu_dec%m)]
    call symblu_sprowcold(temp_lu_dec, piv_order)
  
  else
    write(*,fmt_output) 'Wrong input --> STOP!'
    stop
  end if

  ! Convert back to CSR format
  lu_miter = rowcold_to_csr(temp_lu_dec)

  ! Get the permutation vector and map values of miter to the permuted LU matrix
  call get_lu_permutation(lu_perm, lu_miter, miter)
  
  call cpu_time(time_symbolic)
  write(*,*)
  write(*,'(10X,A,1X,F12.8,A)') 'Time symbolic factorization = ', &
                                 time_symbolic - start_timer, ' [sec]'

  ! Generate random test data
  allocate(rhs(lu_miter%m))
  call random_number(lu_miter%val)
  call random_number(rhs)

  ! Numerical factorization
  call cpu_time(start_timer)
  call sparse_lu(lu_miter)
  call cpu_time(time_numeric)
  write(*,'(10X,A,1X,F12.8,A)') 'Time numeric factorization  = ', &
                                 time_numeric - start_timer, ' [sec]'

  ! Solve the triangular systems
  call cpu_time(start_timer)
  call solve_sparse(lu_miter, rhs)
  call cpu_time(time_solve)
  write(*,'(10X,A,1X,F12.8,A)') 'Time solving triangular sys = ', &
                                 time_solve - start_timer, ' [sec]'

end program main_sparse_gauss_elimination
