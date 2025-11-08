!===============================================================================
! Module: sparse_mod
!
! Purpose: Collection of sparse matrix calculations for chemical reaction systems
!
! Description:
!   This module implements sparse matrix operations using Compressed Sparse Row
!   (CSR) format as the main storage format. Includes symbolic factorization,
!   LU decomposition, matrix operations, and linear system solvers.
!
! Main Features:
!   - Symbolic and numeric matrix-matrix multiplication
!   - Symbolic and numeric matrix addition/subtraction
!   - Matrix transposition
!   - Sparse LU factorization with Markowitz ordering
!   - Forward/backward substitution for triangular systems
!   - Format conversions between different sparse formats
!   - Matrix I/O operations
!
! Author: Willi Schimmel
! Institute: Leibniz Institute for Tropospheric Research (TROPOS)
!===============================================================================
module sparse_mod

  use kind_mod
  use mo_unirnk

  implicit none

  ! Named constants
  integer, parameter, private :: initial_len = 400
  integer, parameter, private :: add_len = 10
  integer, parameter, private :: undefined_ptr = -42
  integer, parameter, private :: uninitialized = -99
  real(dp), parameter, private :: undefined_val = -99999999999999.0_dp
  integer, parameter, private :: default_io_unit = 99

  ! Compressed Sparse Row (CSR) matrix type
  type :: csr_matrix_t
    integer :: m = 0, n = 0, nnz = 0
    integer, allocatable :: row_ptr(:)
    integer, allocatable :: col_ind(:)
    integer, allocatable :: diag_ptr(:)
    integer, allocatable :: diag_ptr_p(:)
    integer, allocatable :: diag_ptr_r(:)
    integer, allocatable :: diag_ptr_c(:)
    integer, allocatable :: row_vector_ptr(:)
    integer, allocatable :: col_vector_ptr(:)
    integer :: x_ptr = undefined_ptr
    integer, allocatable :: permu(:)
    integer, allocatable :: inv_per(:)
    integer, allocatable :: lu_perm(:)
    real(dp), allocatable :: val(:)
  end type csr_matrix_t

  ! Standard row index, column index format
  type :: sparse_row_ind_col_ind_t
    integer :: m = 0, n = 0, nnz = 0
    integer, allocatable :: row_ind(:)
    integer, allocatable :: col_ind(:)
    real(dp), allocatable :: val(:)
  end type sparse_row_ind_col_ind_t

  ! Dynamic sparse row-column format for symbolic factorization
  type :: sprowcold_t
    integer :: m = 0, n = 0
    integer, pointer :: row_ptr(:,:) => null()
    integer, pointer :: col_ind(:) => null()
    integer, pointer :: permu(:) => null()
    integer, pointer :: inv_per(:) => null()
    integer :: ep = 1
    integer :: last = 0
    integer :: len = 0
    integer :: nnz = 0
  end type sprowcold_t

  ! Operator overloading
  interface operator(*)
    module procedure dax_sparse
  end interface

contains

  !-----------------------------------------------------------------------------
  ! Function: new_csr
  !
  ! Purpose: Create and initialize a new CSR matrix
  !
  ! Arguments:
  !   m, n    - matrix dimensions
  !   nnz     - (optional) number of nonzeros
  !   ri, ci  - (optional) row and column indices
  !   val     - (optional) values
  !-----------------------------------------------------------------------------
  function new_csr(m, n, nnz, ri, ci, val) result(new_mat)
    integer, intent(in) :: m, n
    integer, optional, intent(in) :: nnz
    integer, optional, intent(in) :: ri(:), ci(:)
    real(dp), optional, intent(in) :: val(:)
    type(csr_matrix_t) :: new_mat

    integer :: i, same_row, c_cnt

    new_mat%m = m
    new_mat%n = n

    allocate(new_mat%row_ptr(m+1))
    new_mat%row_ptr = 0
    new_mat%row_ptr(1) = 1

    if (present(nnz)) then
      allocate(new_mat%col_ind(nnz))
      new_mat%col_ind = -1
      allocate(new_mat%val(nnz))
      new_mat%val = 0.0_dp
      new_mat%nnz = nnz
    end if

    ! Build CSR structure from row and column indices if provided
    if (present(ri) .and. present(ci)) then
      if (size(ri) /= size(ci)) then
        write(*,'(A)') 'ERROR: Row and column index arrays must have same size'
        stop
      end if
      
      do i = 1, m
        same_row = count(ri == i)
        new_mat%row_ptr(i+1) = new_mat%row_ptr(i) + same_row
      end do
      new_mat%col_ind = ci
      if (present(val)) new_mat%val = val
    end if
  end function new_csr

  !-----------------------------------------------------------------------------
  ! Subroutine: free_matrix_csr
  !
  ! Purpose: Deallocate CSR matrix arrays
  !-----------------------------------------------------------------------------
  subroutine free_matrix_csr(a)
    type(csr_matrix_t), intent(inout) :: a

    if (allocated(a%row_ptr)) deallocate(a%row_ptr)
    if (allocated(a%col_ind)) deallocate(a%col_ind)
    if (allocated(a%diag_ptr)) deallocate(a%diag_ptr)
    if (allocated(a%diag_ptr_r)) deallocate(a%diag_ptr_r)
    if (allocated(a%diag_ptr_c)) deallocate(a%diag_ptr_c)
    if (allocated(a%permu)) deallocate(a%permu)
    if (allocated(a%inv_per)) deallocate(a%inv_per)
    if (allocated(a%val)) deallocate(a%val)
  end subroutine free_matrix_csr

  !-----------------------------------------------------------------------------
  ! Subroutine: free_sprowcold
  !
  ! Purpose: Deallocate dynamic sparse row-column structure
  !-----------------------------------------------------------------------------
  subroutine free_sprowcold(a)
    type(sprowcold_t), intent(inout) :: a

    a%m = 0
    a%n = 0
    a%ep = 1
    a%last = 0
    a%len = 0
    a%nnz = 0
    if (associated(a%row_ptr)) nullify(a%row_ptr)
    if (associated(a%col_ind)) nullify(a%col_ind)
    if (associated(a%permu)) nullify(a%permu)
    if (associated(a%inv_per)) nullify(a%inv_per)
  end subroutine free_sprowcold

  !-----------------------------------------------------------------------------
  ! Subroutine: free_sparse_row_ind_col_ind
  !
  ! Purpose: Deallocate sparse row-column index structure
  !-----------------------------------------------------------------------------
  subroutine free_sparse_row_ind_col_ind(a)
    type(sparse_row_ind_col_ind_t), intent(inout) :: a

    if (allocated(a%row_ind)) deallocate(a%row_ind)
    if (allocated(a%col_ind)) deallocate(a%col_ind)
    if (allocated(a%val)) deallocate(a%val)
  end subroutine free_sparse_row_ind_col_ind

  !-----------------------------------------------------------------------------
  ! Function: sparse_identity
  !
  ! Purpose: Create a sparse identity matrix
  !-----------------------------------------------------------------------------
  function sparse_identity(dim) result(mat)
    integer, intent(in) :: dim
    type(csr_matrix_t) :: mat
    integer :: i

    mat = new_csr(dim, dim, dim)
    do i = 1, dim
      mat%row_ptr(i+1) = mat%row_ptr(i) + 1
      mat%col_ind(i) = i
    end do
    mat%val = 1.0_dp
  end function sparse_identity

  !-----------------------------------------------------------------------------
  ! Function: csr_to_full
  !
  ! Purpose: Convert CSR matrix to full dense format
  !-----------------------------------------------------------------------------
  function csr_to_full(csr) result(full)
    type(csr_matrix_t), intent(in) :: csr
    real(dp) :: full(csr%m, csr%n)
    integer :: i, j, jj

    full = 0.0_dp
    do i = 1, csr%m
      do jj = csr%row_ptr(i), csr%row_ptr(i+1) - 1
        j = csr%col_ind(jj)
        full(i,j) = csr%val(jj)
      end do
    end do
  end function csr_to_full

  !-----------------------------------------------------------------------------
  ! Function: full_to_csr
  !
  ! Purpose: Convert full dense matrix to CSR format
  !-----------------------------------------------------------------------------
  function full_to_csr(full) result(csr)
    real(dp), intent(in) :: full(:,:)
    type(csr_matrix_t) :: csr
    integer :: i, j, jj
    integer :: m, n, nnz

    m = size(full, 1)
    n = size(full, 2)
    nnz = count(full /= 0.0_dp)

    csr = new_csr(m, n, nnz)

    jj = 0
    do i = 1, m
      csr%row_ptr(i+1) = csr%row_ptr(i) + count(full(i,:) /= 0.0_dp)
      do j = 1, n
        if (full(i,j) /= 0.0_dp) then
          jj = jj + 1
          csr%col_ind(jj) = j
          csr%val(jj) = full(i,j)
        end if
      end do
    end do
    csr%nnz = jj
  end function full_to_csr

  !-----------------------------------------------------------------------------
  ! Function: copy_csr
  !
  ! Purpose: Create a deep copy of a CSR matrix
  !-----------------------------------------------------------------------------
  function copy_csr(orig) result(copy)
    type(csr_matrix_t), intent(in) :: orig
    type(csr_matrix_t) :: copy

    copy = new_csr(orig%m, orig%n, orig%nnz)
    copy%row_ptr = orig%row_ptr
    copy%col_ind = orig%col_ind
    copy%val = orig%val
    copy%nnz = orig%nnz
    
    if (orig%x_ptr /= undefined_ptr) copy%x_ptr = orig%x_ptr
    if (allocated(orig%diag_ptr)) copy%diag_ptr = orig%diag_ptr
    if (allocated(orig%diag_ptr_p)) copy%diag_ptr_p = orig%diag_ptr_p
    if (allocated(orig%diag_ptr_r)) copy%diag_ptr_r = orig%diag_ptr_r
    if (allocated(orig%diag_ptr_c)) copy%diag_ptr_c = orig%diag_ptr_c
    if (allocated(orig%row_vector_ptr)) copy%row_vector_ptr = orig%row_vector_ptr
    if (allocated(orig%col_vector_ptr)) copy%col_vector_ptr = orig%col_vector_ptr
    if (allocated(orig%permu)) copy%permu = orig%permu
    if (allocated(orig%inv_per)) copy%inv_per = orig%inv_per
    if (allocated(orig%lu_perm)) copy%lu_perm = orig%lu_perm
  end function copy_csr

  !-----------------------------------------------------------------------------
  ! Function: rowcold_to_csr
  !
  ! Purpose: Convert dynamic row-column format to CSR format
  !-----------------------------------------------------------------------------
  function rowcold_to_csr(sp_row) result(csr)
    type(sprowcold_t), intent(in) :: sp_row
    type(csr_matrix_t) :: csr
    integer :: i, j, jj, nzr_a

    csr = new_csr(sp_row%m, sp_row%n, sp_row%nnz)
    allocate(csr%diag_ptr(sp_row%m), csr%diag_ptr_p(sp_row%m))
    csr%diag_ptr = -11
    csr%diag_ptr_p = -11

    nzr_a = 0
    do i = 1, csr%m
      csr%row_ptr(i+1) = csr%row_ptr(i)
      do jj = sp_row%row_ptr(1,i), sp_row%row_ptr(2,i)
        j = sp_row%col_ind(jj)
        nzr_a = nzr_a + 1
        csr%col_ind(nzr_a) = j
        if (i == j) then
          csr%diag_ptr(i) = csr%row_ptr(i+1)
        end if
        csr%row_ptr(i+1) = csr%row_ptr(i+1) + 1
      end do
    end do

    if (associated(sp_row%permu)) then
      if (.not. allocated(csr%permu)) then
        allocate(csr%permu(csr%n))
        csr%permu(:) = -16
      end if
      csr%permu(:) = sp_row%permu(:)
    end if
    
    if (associated(sp_row%inv_per)) then
      if (.not. allocated(csr%inv_per)) then
        allocate(csr%inv_per(csr%n))
        csr%inv_per(:) = -16
      end if
      csr%inv_per(:) = sp_row%inv_per(:)
    end if
  end function rowcold_to_csr

  !-----------------------------------------------------------------------------
  ! Function: csr_to_sprowcold
  !
  ! Purpose: Convert CSR to dynamic row-column format for factorization
  !-----------------------------------------------------------------------------
  function csr_to_sprowcold(csr) result(sp_row_col)
    type(csr_matrix_t), intent(in) :: csr
    type(sprowcold_t) :: sp_row_col
    integer :: nzr_csr, start_idx, end_idx
    integer :: i, j, jj

    sp_row_col%m = csr%m
    sp_row_col%n = csr%n
    allocate(sp_row_col%row_ptr(2, sp_row_col%n))
    
    nzr_csr = size(csr%col_ind)
    sp_row_col%len = 10 * nzr_csr + add_len * sp_row_col%n
    allocate(sp_row_col%col_ind(sp_row_col%len))
    sp_row_col%col_ind = 0
    allocate(sp_row_col%permu(sp_row_col%n))
    sp_row_col%permu = 0
    allocate(sp_row_col%inv_per(sp_row_col%n))
    sp_row_col%inv_per = 0

    start_idx = 1
    do i = 1, csr%n
      sp_row_col%row_ptr(1,i) = start_idx
      end_idx = start_idx + (csr%row_ptr(i+1) - csr%row_ptr(i) - 1)
      sp_row_col%row_ptr(2,i) = end_idx
      
      do jj = csr%row_ptr(i), csr%row_ptr(i+1) - 1
        j = csr%col_ind(jj)
        sp_row_col%col_ind(start_idx) = j
        start_idx = start_idx + 1
      end do
      start_idx = start_idx + add_len - 1
    end do
    
    sp_row_col%ep = sp_row_col%row_ptr(2, sp_row_col%n) + 1
    sp_row_col%last = sp_row_col%n
  end function csr_to_sprowcold

  !-----------------------------------------------------------------------------
  ! Subroutine: sort_vec
  !
  ! Purpose: Simple bubble sort for integer vectors (small vectors only)
  !-----------------------------------------------------------------------------
  subroutine sort_vec(vec)
    integer, intent(inout) :: vec(:)
    integer :: i, j, n, i_temp

    n = size(vec)
    do i = 1, n
      do j = 1, n - i
        if (vec(j) > vec(j+1)) then
          i_temp = vec(j)
          vec(j) = vec(j+1)
          vec(j+1) = i_temp
        end if
      end do
    end do
  end subroutine sort_vec

  !-----------------------------------------------------------------------------
  ! Subroutine: swap_int
  !
  ! Purpose: Swap two integer values
  !-----------------------------------------------------------------------------
  subroutine swap_int(i, j)
    integer, intent(inout) :: i, j
    integer :: i_temp
    
    i_temp = i
    i = j
    j = i_temp
  end subroutine swap_int

  !-----------------------------------------------------------------------------
  ! Subroutine: symblu_sprowcold
  !
  ! Purpose: Symbolic LU factorization with given pivot ordering
  !-----------------------------------------------------------------------------
  subroutine symblu_sprowcold(a, permu)
    type(sprowcold_t), intent(inout) :: a
    integer, intent(in) :: permu(:)
    
    integer :: row_piv(a%n)
    integer :: i, j, l, jj, ip, i_piv
    logical :: ins

    do i = 1, a%n
      a%inv_per(i) = i
      a%permu(i) = i
    end do

    do i = 1, a%n
      ip = a%permu(permu(i))
      call swap_int(a%inv_per(i), a%inv_per(ip))
      call swap_int(a%row_ptr(1,i), a%row_ptr(1,ip))
      call swap_int(a%row_ptr(2,i), a%row_ptr(2,ip))
      a%permu(a%inv_per(i)) = i
      a%permu(a%inv_per(ip)) = ip
      
      if (a%last == i) then
        a%last = ip
      else if (a%last == ip) then
        a%last = i
      end if

      ! Update: find pivot row elements
      i_piv = 0
      do jj = a%row_ptr(1,i), a%row_ptr(2,i)
        if (a%permu(a%col_ind(jj)) > i) then
          i_piv = i_piv + 1
          row_piv(i_piv) = a%col_ind(jj)
        end if
      end do
      
      ! Insert fill-in elements
      if (i_piv > 0) then
        do j = i+1, a%n
          do jj = a%row_ptr(1,j), a%row_ptr(2,j)
            if (a%permu(a%col_ind(jj)) == i) then
              do l = 1, i_piv
                call insert_sprowcold(a, j, row_piv(l), ins)
              end do
              exit
            end if
          end do
        end do
      end if
    end do
    
    ! Apply permutation and sort
    do i = 1, a%n
      do jj = a%row_ptr(1,i), a%row_ptr(2,i)
        a%col_ind(jj) = a%permu(a%col_ind(jj))
      end do
      call sort_vec(a%col_ind(a%row_ptr(1,i):a%row_ptr(2,i)))
      a%nnz = a%nnz + size(a%col_ind(a%row_ptr(1,i):a%row_ptr(2,i)))
    end do
  end subroutine symblu_sprowcold

  !-----------------------------------------------------------------------------
  ! Subroutine: symblu_sprowcold_m
  !
  ! Purpose: Symbolic LU factorization with Markowitz minimum-degree ordering
  !-----------------------------------------------------------------------------
  subroutine symblu_sprowcold_m(a)
    type(sprowcold_t), intent(inout) :: a
    
    integer :: r(a%n), c(a%n), row_piv(a%n)
    integer :: i, j, l, jj, ip, i_piv
    real(dp) :: md
    integer :: rc
    logical :: ins

    c = 0
    do i = 1, a%n
      a%inv_per(i) = i
      a%permu(i) = i
      r(i) = a%row_ptr(2,i) - a%row_ptr(1,i) + 1
      do jj = a%row_ptr(1,i), a%row_ptr(2,i)
        c(a%col_ind(jj)) = c(a%col_ind(jj)) + 1
      end do
    end do

    ! Main Markowitz ordering loop
    do i = 1, a%n
      ip = 0
      md = 1.0e99_dp

      ! Find pivot with minimum Markowitz cost
      do j = i, a%n
        rc = (r(j) - 1) * (c(j) - 1)
        if (rc <= md) then
          md = real(rc, dp)
          ip = j
        end if
      end do

      call swap_int(r(i), r(ip))
      call swap_int(c(i), c(ip))
      call swap_int(a%inv_per(i), a%inv_per(ip))
      call swap_int(a%row_ptr(1,i), a%row_ptr(1,ip))
      call swap_int(a%row_ptr(2,i), a%row_ptr(2,ip))

      a%permu(a%inv_per(i)) = i
      a%permu(a%inv_per(ip)) = ip

      if (a%last == i) then
        a%last = ip
      else if (a%last == ip) then
        a%last = i
      end if

      ! Update Markowitz counts and structure
      i_piv = 0
      do jj = a%row_ptr(1,i), a%row_ptr(2,i)
        if (a%permu(a%col_ind(jj)) > i) then
          i_piv = i_piv + 1
          row_piv(i_piv) = a%col_ind(jj)
          c(a%permu(a%col_ind(jj))) = c(a%permu(a%col_ind(jj))) - 1
        end if
      end do
      
      if (i_piv > 0) then
        do j = i+1, a%n
          do jj = a%row_ptr(1,j), a%row_ptr(2,j)
            if (a%permu(a%col_ind(jj)) == i) then
              r(j) = r(j) - 1
              c(i) = c(i) - 1
              do l = 1, i_piv
                call insert_sprowcold(a, j, row_piv(l), ins)
                if (ins) then
                  c(a%permu(row_piv(l))) = c(a%permu(row_piv(l))) + 1
                  r(j) = r(j) + 1
                end if
              end do
              exit
            end if
          end do
          if (c(i) == 1) exit
        end do
      end if
    end do

    ! Apply permutation and sort
    do i = 1, a%n
      do jj = a%row_ptr(1,i), a%row_ptr(2,i)
        a%col_ind(jj) = a%permu(a%col_ind(jj))
      end do
      call sort_vec(a%col_ind(a%row_ptr(1,i):a%row_ptr(2,i)))
      a%nnz = a%nnz + size(a%col_ind(a%row_ptr(1,i):a%row_ptr(2,i)))
    end do
  end subroutine symblu_sprowcold_m

  !-----------------------------------------------------------------------------
  ! Subroutine: insert_sprowcold
  !
  ! Purpose: Insert element into dynamic sparse structure
  !-----------------------------------------------------------------------------
  subroutine insert_sprowcold(a, i_a, j_a, ins)
    type(sprowcold_t), intent(inout) :: a
    integer, intent(in) :: i_a, j_a
    logical, intent(out) :: ins
    
    integer :: i_temp, j, l
    integer, allocatable :: i_work(:)

    if (i_a == 0 .or. j_a == 0) then
      write(*,'(A,I0)') 'ERROR in insert_sprowcold: i_a = ', i_a
      write(*,'(A,I0)') 'ERROR in insert_sprowcold: j_a = ', j_a
      stop
    end if

    ! Check if element already exists
    ins = .true.
    do j = a%row_ptr(1,i_a), a%row_ptr(2,i_a)
      if (j_a == a%col_ind(j)) then
        ins = .false.
        return
      end if
    end do

    ! Insert new element
    if (ins) then
      if (a%col_ind(a%row_ptr(2,i_a)+1) /= 0) then
        ! Move row to end of storage
        i_temp = a%ep
        do l = a%row_ptr(1,i_a), a%row_ptr(2,i_a)
          a%col_ind(a%ep) = a%col_ind(l)
          a%col_ind(l) = 0
          a%ep = a%ep + 1
        end do
        a%row_ptr(2,i_a) = a%ep - 1
        a%row_ptr(1,i_a) = i_temp
        a%last = i_a
      end if
      
      a%row_ptr(2,i_a) = a%row_ptr(2,i_a) + 1
      a%col_ind(a%row_ptr(2,i_a)) = j_a
      if (i_a == a%last) a%ep = a%ep + 1
    end if

    ! Garbage collection if needed
    if (a%ep >= a%len - a%m) then
      call gcmat_sprowcold(a)
    end if

    ! Expand storage if needed
    if (a%ep >= a%len - a%m) then
      allocate(i_work(a%ep))
      i_work(1:a%ep) = a%col_ind(1:a%ep)
      deallocate(a%col_ind)
      a%len = 2 * a%len
      allocate(a%col_ind(a%len))
      a%col_ind(1:a%ep) = i_work(1:a%ep)
      deallocate(i_work)
    end if
  end subroutine insert_sprowcold

  !-----------------------------------------------------------------------------
  ! Subroutine: gcmat_sprowcold
  !
  ! Purpose: Garbage collection for dynamic sparse structure
  !           Compacts the storage by removing holes
  !-----------------------------------------------------------------------------
  subroutine gcmat_sprowcold(a)
    type(sprowcold_t), intent(inout) :: a
    
    integer :: i, iz, j, l, pointr, rowlen, ep, len, m

    m = a%m
    ep = a%ep
    len = a%len
    pointr = 1
    i = 1

    do
      if (i >= ep) exit
      
      if (a%col_ind(i) /= 0) then
        ! Find current row and its length
        do l = 1, m
          if (a%row_ptr(1,l) <= i .and. i <= a%row_ptr(2,l)) then
            iz = l
            exit
          end if
        end do
        rowlen = a%row_ptr(2,iz) - a%row_ptr(1,iz)

        ! Set new addresses for current row
        a%row_ptr(1,iz) = pointr
        a%row_ptr(2,iz) = pointr + rowlen
        
        do j = pointr, pointr + rowlen
          a%col_ind(j) = a%col_ind(i)
          i = i + 1
        end do
        i = i - 1
        pointr = a%row_ptr(2,iz) + 1
      end if
      i = i + 1
    end do

    ! Set free space
    ep = pointr
    do i = 1, m
      if (a%row_ptr(1,i) > ep) then
        a%row_ptr(1,i) = ep
        a%row_ptr(2,i) = a%row_ptr(1,i) - 1
        ep = ep + initial_len
      end if
    end do

    do i = pointr, len
      a%col_ind(i) = 0
    end do
    a%ep = ep
  end subroutine gcmat_sprowcold

  !-----------------------------------------------------------------------------
  ! Subroutine: transpose_sparse
  !
  ! Purpose: Compute transpose of a sparse matrix
  !-----------------------------------------------------------------------------
  subroutine transpose_sparse(mat_at, mat_a)
    type(csr_matrix_t), intent(in) :: mat_a
    type(csr_matrix_t), intent(out) :: mat_at
    integer :: i, j, indx

    mat_at = new_csr(mat_a%n, mat_a%m)

    ! Count entries per row in transpose
    do i = 1, mat_a%m
      do j = mat_a%row_ptr(i), mat_a%row_ptr(i+1) - 1
        mat_at%row_ptr(mat_a%col_ind(j)+1) = mat_at%row_ptr(mat_a%col_ind(j)+1) + 1
      end do
    end do

    ! Cumulative sum to get row pointers
    do i = 1, mat_at%m
      mat_at%row_ptr(i+1) = mat_at%row_ptr(i) + mat_at%row_ptr(i+1)
    end do

    ! Fill in values
    allocate(mat_at%col_ind(size(mat_a%col_ind)))
    allocate(mat_at%val(size(mat_a%val)))
    mat_at%col_ind = 0
    mat_at%val = 0.0_dp
    
    do i = 1, mat_a%m
      do j = mat_a%row_ptr(i), mat_a%row_ptr(i+1) - 1
        indx = mat_a%col_ind(j)
        mat_at%col_ind(mat_at%row_ptr(indx)) = i
        mat_at%val(mat_at%row_ptr(indx)) = mat_a%val(j)
        mat_at%row_ptr(indx) = mat_at%row_ptr(indx) + 1
      end do
    end do
    
    ! Restore row pointers
    do i = mat_at%m, 1, -1
      mat_at%row_ptr(i+1) = mat_at%row_ptr(i)
    end do
    mat_at%row_ptr(1) = 1
    mat_at%nnz = mat_a%row_ptr(mat_a%m+1) - 1
  end subroutine transpose_sparse

  !-----------------------------------------------------------------------------
  ! Subroutine: symbolic_mult
  !
  ! Purpose: Symbolic matrix multiplication C = A * B
  !          Computes the sparsity pattern without computing values
  !-----------------------------------------------------------------------------
  subroutine symbolic_mult(a, b, c)
    type(csr_matrix_t), intent(in) :: a, b
    type(csr_matrix_t), intent(out) :: c
    
    integer :: indx(max(a%n, a%m, b%n))
    integer :: i, j, jj, k
    integer :: i_start, length, i_temp

    c = new_csr(a%m, b%n)

    ! First pass: determine row lengths
    indx = 0
    do i = 1, a%m
      i_start = -1
      length = 0
      do jj = a%row_ptr(i), a%row_ptr(i+1) - 1
        j = a%col_ind(jj)
        do k = b%row_ptr(j), b%row_ptr(j+1) - 1
          if (indx(b%col_ind(k)) == 0) then
            indx(b%col_ind(k)) = i_start
            i_start = b%col_ind(k)
            length = length + 1
          end if
        end do
      end do
      c%row_ptr(i+1) = c%row_ptr(i) + length

      ! Clear index array
      do j = c%row_ptr(i), c%row_ptr(i+1) - 1
        i_temp = i_start
        i_start = indx(i_start)
        indx(i_temp) = 0
      end do
      indx(i) = 0
    end do

    ! Allocate and fill column indices
    allocate(c%col_ind(c%row_ptr(c%m+1) - 1))
    c%col_ind = 0

    indx = 0
    do i = 1, a%m
      i_start = -1
      length = 0
      do jj = a%row_ptr(i), a%row_ptr(i+1) - 1
        j = a%col_ind(jj)
        do k = b%row_ptr(j), b%row_ptr(j+1) - 1
          if (indx(b%col_ind(k)) == 0) then
            indx(b%col_ind(k)) = i_start
            i_start = b%col_ind(k)
            length = length + 1
          end if
        end do
      end do
      c%row_ptr(i+1) = c%row_ptr(i) + length
      
      do j = c%row_ptr(i), c%row_ptr(i+1) - 1
        c%col_ind(j) = i_start
        i_start = indx(i_start)
        indx(c%col_ind(j)) = 0
      end do
      indx(i) = 0
    end do

    ! Sort column indices
    do i = 1, c%m
      call sort_vec(c%col_ind(c%row_ptr(i):c%row_ptr(i+1)-1))
    end do
    
    allocate(c%val(c%row_ptr(c%m+1) - 1))
    c%val = 0.0_dp
    c%nnz = c%row_ptr(c%m+1) - 1
  end subroutine symbolic_mult

  !-----------------------------------------------------------------------------
  ! Subroutine: sparse_mult
  !
  ! Purpose: Numeric matrix multiplication C = A * B (structure already known)
  !-----------------------------------------------------------------------------
  subroutine sparse_mult(a, b, c)
    type(csr_matrix_t), intent(inout) :: a, c
    type(csr_matrix_t), intent(in) :: b
    
    integer :: i, j, jj, k, kk
    real(dp) :: ajj
    real(dp) :: temp(max(a%m, a%n, b%n))

    temp = 0.0_dp

    do i = 1, a%m
      do jj = a%row_ptr(i), a%row_ptr(i+1) - 1
        j = a%col_ind(jj)
        ajj = a%val(jj)
        do kk = b%row_ptr(j), b%row_ptr(j+1) - 1
          k = b%col_ind(kk)
          temp(k) = temp(k) + ajj * b%val(kk)
        end do
      end do
      
      do j = c%row_ptr(i), c%row_ptr(i+1) - 1
        c%val(j) = temp(c%col_ind(j))
        temp(c%col_ind(j)) = 0.0_dp
      end do
    end do
  end subroutine sparse_mult

  !-----------------------------------------------------------------------------
  ! Subroutine: symbolic_add
  !
  ! Purpose: Symbolic matrix addition C = A + B
  !          Computes sparsity pattern without computing values
  !-----------------------------------------------------------------------------
  subroutine symbolic_add(mat_c, mat_a, mat_b)
    type(csr_matrix_t), intent(in) :: mat_a, mat_b
    type(csr_matrix_t), intent(out) :: mat_c
    
    integer :: i, ii, j, jj, k, kk
    integer, allocatable :: tmp_col(:)
    integer, allocatable :: perm_vec(:)
    integer :: col_len, current_length, same_cnt

    if (.not.((mat_a%m == mat_b%m) .and. (mat_a%n == mat_b%n))) then
      write(*,'(A)') 'ERROR: Matrix dimensions do not match'
      write(*,'(A,2I6)') 'Matrix A: ', mat_a%m, mat_a%n
      write(*,'(A,2I6)') 'Matrix B: ', mat_b%m, mat_b%n
      stop
    end if

    mat_c = new_csr(mat_a%m, mat_a%n)

    ! Determine row lengths
    do i = 1, mat_c%m
      current_length = (mat_a%row_ptr(i+1) - mat_a%row_ptr(i)) + &
                      (mat_b%row_ptr(i+1) - mat_b%row_ptr(i))
      same_cnt = 0
      
      do ii = mat_a%row_ptr(i), mat_a%row_ptr(i+1) - 1
        do j = mat_b%row_ptr(i), mat_b%row_ptr(i+1) - 1
          if (mat_a%col_ind(ii) == mat_b%col_ind(j)) then
            same_cnt = same_cnt + 1
          end if
        end do
      end do
      mat_c%row_ptr(i+1) = mat_c%row_ptr(i) + current_length - same_cnt
    end do

    ! Allocate and fill column indices
    allocate(mat_c%col_ind(mat_c%row_ptr(mat_c%m+1) - 1))
    mat_c%col_ind = 0

    kk = 1
    do i = 1, mat_c%m
      k = 1
      current_length = mat_c%row_ptr(i+1) - mat_c%row_ptr(i)
      same_cnt = 0
      
      do ii = mat_a%row_ptr(i), mat_a%row_ptr(i+1) - 1
        do j = mat_b%row_ptr(i), mat_b%row_ptr(i+1) - 1
          if (mat_a%col_ind(ii) == mat_b%col_ind(j)) then
            same_cnt = same_cnt + 1
          end if
        end do
      end do

      allocate(tmp_col(current_length + same_cnt))
      tmp_col = 0
      allocate(perm_vec(current_length + same_cnt))
      perm_vec = 0

      do jj = mat_a%row_ptr(i), mat_a%row_ptr(i+1) - 1
        tmp_col(k) = mat_a%col_ind(jj)
        k = k + 1
      end do
      
      do jj = mat_b%row_ptr(i), mat_b%row_ptr(i+1) - 1
        tmp_col(k) = mat_b%col_ind(jj)
        k = k + 1
      end do

      call unirnk(tmp_col, perm_vec, col_len)
      do j = 1, col_len
        mat_c%col_ind(kk) = tmp_col(perm_vec(j))
        kk = kk + 1
      end do

      deallocate(tmp_col)
      deallocate(perm_vec)
    end do

    allocate(mat_c%val(mat_c%row_ptr(mat_c%m+1) - 1))
    mat_c%val = 0.0_dp
    mat_c%nnz = mat_c%row_ptr(mat_c%m+1) - 1
  end subroutine symbolic_add

  !-----------------------------------------------------------------------------
  ! Subroutine: sparse_add
  !
  ! Purpose: Numeric matrix addition/subtraction C = A +/- B
  !-----------------------------------------------------------------------------
  subroutine sparse_add(mat_c, mat_a, mat_b, sub)
    type(csr_matrix_t), intent(in) :: mat_a, mat_b
    type(csr_matrix_t), intent(inout) :: mat_c
    character, optional, intent(in) :: sub
    
    real(dp), allocatable :: tmp_b_val(:)
    integer :: i, ii, jj, kk

    allocate(tmp_b_val(mat_b%row_ptr(mat_b%m+1) - 1))
    if (present(sub)) then
      if (sub == '-') then
        tmp_b_val = -mat_b%val
      else
        tmp_b_val = mat_b%val
      end if
    else
      tmp_b_val = mat_b%val
    end if

    do i = 1, mat_c%m
      do ii = mat_c%row_ptr(i), mat_c%row_ptr(i+1) - 1
        do jj = mat_a%row_ptr(i), mat_a%row_ptr(i+1) - 1
          if (mat_a%col_ind(jj) == mat_c%col_ind(ii)) then
            mat_c%val(ii) = mat_c%val(ii) + mat_a%val(jj)
          end if
        end do
        do kk = mat_b%row_ptr(i), mat_b%row_ptr(i+1) - 1
          if (mat_b%col_ind(kk) == mat_c%col_ind(ii)) then
            mat_c%val(ii) = mat_c%val(ii) + tmp_b_val(kk)
          end if
        end do
      end do
    end do
    deallocate(tmp_b_val)
  end subroutine sparse_add

  !-----------------------------------------------------------------------------
  ! Subroutine: get_lu_permutation
  !
  ! Purpose: Get permutation vector mapping original to LU structure
  !-----------------------------------------------------------------------------
  subroutine get_lu_permutation(permutation, lu, a)
    type(csr_matrix_t), intent(inout) :: lu, a
    integer, allocatable, intent(out) :: permutation(:)
    
    integer :: nnz_a
    integer :: i, ip, j, jj, jp, jjp, jp1

    nnz_a = a%row_ptr(a%m+1) - 1

    if (.not. allocated(permutation)) allocate(permutation(nnz_a))
    permutation = -14

    lu%val = 0.0_dp

    do i = 1, a%n
      ip = lu%permu(i)
      do jj = a%row_ptr(i), a%row_ptr(i+1) - 1
        j = a%col_ind(jj)
        jp = lu%permu(a%col_ind(jj))
        do jjp = lu%row_ptr(ip), lu%row_ptr(ip+1) - 1
          jp1 = lu%col_ind(jjp)
          if (jp1 == jp) then
            lu%val(jjp) = a%val(jj)
            permutation(jj) = jjp
          end if
        end do
      end do
    end do
    
    allocate(a%lu_perm(nnz_a), lu%lu_perm(nnz_a))
    a%lu_perm = permutation
    lu%lu_perm = permutation
  end subroutine get_lu_permutation

  !-----------------------------------------------------------------------------
  ! Subroutine: print_sparse_matrix
  !
  ! Purpose: Print sparse matrix to stdout
  !-----------------------------------------------------------------------------
  subroutine print_sparse_matrix(a)
    type(csr_matrix_t), intent(in) :: a
    integer :: i, j, jj

    write(*,'(A,2I6)') 'Dimension:  ', a%m, a%n
    write(*,'(A,I8)') 'Nonzeros:   ', size(a%col_ind)
    write(*,*)

    do i = 1, a%m
      do jj = a%row_ptr(i), a%row_ptr(i+1) - 1
        j = a%col_ind(jj)
        write(*,'(1X,I5,1X,I5,10X,ES23.14)') i, j, a%val(jj)
      end do
    end do
    write(*,*)
  end subroutine print_sparse_matrix

  !-----------------------------------------------------------------------------
  ! Subroutine: write_sparse_matrix
  !
  ! Purpose: Write sparse matrix to file
  !-----------------------------------------------------------------------------
  subroutine write_sparse_matrix(a, filename, nr, ns, extended, classic, teq, linalg)
    type(csr_matrix_t), intent(in) :: a
    character(*), optional, intent(in) :: filename
    integer, optional, intent(in) :: nr, ns
    logical, optional, intent(in) :: extended, classic, teq
    character(2), optional, intent(in) :: linalg
    
    integer :: i, j, jj
    character(len=256) :: full_filename

    if (present(filename) .and. present(linalg)) then
      full_filename = trim(adjustl(filename)) // '_' // trim(linalg) // '.SparseMat'
    else if (present(filename)) then
      full_filename = trim(adjustl(filename)) // '.SparseMat'
    else
      full_filename = 'matrix.SparseMat'
    end if

    open(unit=default_io_unit, file=trim(full_filename), status='unknown')

    write(default_io_unit,'(A)') '# Sparse Matrix File'
    write(default_io_unit,'(A)') '#'
    if (present(filename)) write(default_io_unit,'(A,A)') '# Name:      ', trim(filename)
    write(default_io_unit,'(A,I0,A,I0)') '# Dimension: ', a%m, ' x ', a%n
    write(default_io_unit,'(A,I0)') '# Nonzeros:  ', a%nnz
    if (present(nr) .and. present(ns)) then
      write(default_io_unit,'(A,I0,A,I0)') '# nreac, nspc: ', nr, ' , ', ns
    end if
    if (present(linalg)) write(default_io_unit,'(A,A)') '# Matrix Form: ', linalg
    write(default_io_unit,'(A)') '#'
    write(default_io_unit,*)
    write(default_io_unit,'(A)') 'MATRIX'
    
    do i = 1, a%m
      do jj = a%row_ptr(i), a%row_ptr(i+1) - 1
        j = a%col_ind(jj)
        write(default_io_unit,'(1X,I12,10X,I12,10X,ES23.14)') i, j, a%val(jj)
      end do
    end do
    write(default_io_unit,*)

    if (allocated(a%diag_ptr)) then
      write(default_io_unit,'(A)') 'DIAG_PTR'
      do i = 1, size(a%diag_ptr)
        write(default_io_unit,'(1X,I6,10X,I6)') i, a%diag_ptr(i)
      end do
      write(default_io_unit,*)
    end if

    if (allocated(a%permu)) then
      write(default_io_unit,'(A)') 'DIAG_PERMUTATION'
      do i = 1, size(a%permu)
        write(default_io_unit,'(1X,I6,10X,I6,10X,I6)') i, a%permu(i), a%inv_per(i)
      end do
      write(default_io_unit,*)
    end if

    if (allocated(a%lu_perm)) then
      write(default_io_unit,'(A)') 'LU_PERMUTATION'
      do i = 1, size(a%lu_perm)
        write(default_io_unit,'(1X,I6,10X,I6)') i, a%lu_perm(i)
      end do
      write(default_io_unit,*)
    end if

    close(default_io_unit)
    write(*,'(A,A)') '  Writing matrix to file: ', trim(full_filename)
  end subroutine write_sparse_matrix

  !-----------------------------------------------------------------------------
  ! Function: read_sparse_matrix
  !
  ! Purpose: Read sparse matrix from file
  !
  ! File Format:
  !   Line 1: m n [comment]
  !   Line 2: nnz [comment]
  !   Line 3: blank
  !   Line 4+: i j value
  !-----------------------------------------------------------------------------
  function read_sparse_matrix(filename) result(a)
    character(*), intent(in) :: filename
    type(csr_matrix_t) :: a
    
    character(len=400) :: io_msg, dummy
    integer :: io_stat
    integer :: i, m, n, nnz
    integer, allocatable :: i_a(:), j_a(:)
    real(dp), allocatable :: v_a(:)

    open(unit=default_io_unit, file=trim(adjustl(filename)), status='old', &
         iostat=io_stat, iomsg=io_msg)
    
    if (io_stat /= 0) then
      write(*,'(A,A)') 'ERROR opening file: ', trim(filename)
      write(*,'(A,A)') 'Message: ', trim(io_msg)
      stop
    end if

    read(default_io_unit,*,iostat=io_stat,iomsg=io_msg) m, n, dummy
    read(default_io_unit,*,iostat=io_stat,iomsg=io_msg) nnz, dummy
    read(default_io_unit,*,iostat=io_stat,iomsg=io_msg) dummy

    allocate(i_a(nnz), j_a(nnz), v_a(nnz))
    i_a = uninitialized
    j_a = uninitialized
    v_a = 0.0_dp

    do i = 1, nnz
      read(default_io_unit,*,iostat=io_stat,iomsg=io_msg) i_a(i), j_a(i), v_a(i)
      if (io_stat /= 0) then
        write(*,'(A,I0)') 'ERROR reading matrix entry ', i
        stop
      end if
    end do

    a = new_csr(m, n, nnz, i_a, j_a, v_a)

    close(default_io_unit)
    write(*,*)
    write(*,'(A,A,A)') '  Read matrix  ::  ', trim(filename), '  done.'
  end function read_sparse_matrix

  !-----------------------------------------------------------------------------
  ! Function: input_pivot_order
  !
  ! Purpose: Read pivot ordering from file
  !-----------------------------------------------------------------------------
  function input_pivot_order(filename, dim) result(piv_ord)
    character(*), intent(in) :: filename
    integer, intent(in) :: dim
    integer, allocatable :: piv_ord(:)
    
    character(len=400) :: io_msg
    integer :: io_stat, i

    open(unit=default_io_unit, file=trim(adjustl(filename)), status='old', &
         iostat=io_stat, iomsg=io_msg)
    
    if (io_stat /= 0) then
      write(*,'(A,A)') 'ERROR opening file: ', trim(filename)
      stop
    end if

    allocate(piv_ord(dim))
    piv_ord = uninitialized

    do i = 1, dim
      read(default_io_unit,*,iostat=io_stat,iomsg=io_msg) piv_ord(i)
    end do

    close(default_io_unit)
    write(*,*)
    write(*,'(A,A,A)') '  Read pivot order  ::  ', trim(filename), '  done.'
  end function input_pivot_order

  !-----------------------------------------------------------------------------
  ! Function: dax_sparse
  !
  ! Purpose: Sparse matrix-vector product: y = A * x
  !-----------------------------------------------------------------------------
  pure function dax_sparse(a, x) result(rhs)
    type(csr_matrix_t), intent(in) :: a
    real(dp), intent(in) :: x(:)
    real(dp) :: rhs(a%m)
    integer :: i, rp_i, rp_i1

    rhs = 0.0_dp
    do i = 1, a%m
      rp_i = a%row_ptr(i)
      rp_i1 = a%row_ptr(i+1) - 1
      rhs(i) = sum(a%val(rp_i:rp_i1) * x(a%col_ind(rp_i:rp_i1)))
    end do
  end function dax_sparse

  !-----------------------------------------------------------------------------
  ! Function: daxpy_sparse
  !
  ! Purpose: Sparse matrix-vector product with addition: y = A * x + y
  !-----------------------------------------------------------------------------
  pure function daxpy_sparse(a, x, y) result(rhs)
    type(csr_matrix_t), intent(in) :: a
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: rhs(a%m)
    integer :: i, rp_i, rp_i1

    rhs = 0.0_dp
    do i = 1, a%m
      rp_i = a%row_ptr(i)
      rp_i1 = a%row_ptr(i+1) - 1
      rhs(i) = sum(a%val(rp_i:rp_i1) * x(a%col_ind(rp_i:rp_i1))) + y(i)
    end do
  end function daxpy_sparse

  !-----------------------------------------------------------------------------
  ! Subroutine: sparse_lu
  !
  ! Purpose: Numeric LU factorization (in-place)
  !          Assumes symbolic factorization has been done
  !-----------------------------------------------------------------------------
  pure subroutine sparse_lu(a)
    type(csr_matrix_t), intent(inout) :: a
    
    real(dp) :: w(a%n)
    real(dp) :: alpha
    integer :: i, j, jj, kk

    do i = 1, a%n
      ! Copy row to working array
      do jj = a%row_ptr(i), a%row_ptr(i+1) - 1
        w(a%col_ind(jj)) = a%val(jj)
      end do
      
      ! Eliminate previous rows
      do jj = a%row_ptr(i), a%diag_ptr(i) - 1
        j = a%col_ind(jj)
        alpha = w(j) / a%val(a%diag_ptr(j))
        w(j) = alpha
        do kk = a%diag_ptr(j) + 1, a%row_ptr(j+1) - 1
          w(a%col_ind(kk)) = w(a%col_ind(kk)) - alpha * a%val(kk)
        end do
      end do
      
      ! Write back results
      do jj = a%row_ptr(i), a%row_ptr(i+1) - 1
        a%val(jj) = w(a%col_ind(jj))
      end do
    end do
  end subroutine sparse_lu

  !-----------------------------------------------------------------------------
  ! Subroutine: solve_sparse
  !
  ! Purpose: Solve LU * x = b using forward/backward substitution
  !          b is overwritten with solution x
  !-----------------------------------------------------------------------------
  pure subroutine solve_sparse(lu, rhs)
    type(csr_matrix_t), intent(in) :: lu
    real(dp), intent(inout) :: rhs(:)
    
    integer :: i, jj
    real(dp) :: b(lu%n)

    ! Apply permutation
    b(lu%permu) = rhs

    ! Forward substitution (L solve)
    do i = 2, lu%n
      do jj = lu%row_ptr(i), lu%diag_ptr(i) - 1
        b(i) = b(i) - lu%val(jj) * b(lu%col_ind(jj))
      end do
    end do

    ! Backward substitution (U solve)
    do i = lu%n, 1, -1
      do jj = lu%diag_ptr(i) + 1, lu%row_ptr(i+1) - 1
        b(i) = b(i) - lu%val(jj) * b(lu%col_ind(jj))
      end do
      b(i) = b(i) / lu%val(lu%diag_ptr(i))
    end do

    ! Reverse permutation
    rhs(lu%inv_per) = b
  end subroutine solve_sparse

end module sparse_mod
