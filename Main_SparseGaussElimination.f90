


PROGRAM Main_SparseGaussElimmination
  !
  USE Kind_Mod
  USE Sparse_Mod

  IMPLICIT NONE
  !
  CHARACTER(80)   :: Filename0 = ''        ! *.run file
  CHARACTER(80)   :: input = ''        ! *.run file


  ! matricies for symbolic phase
  TYPE(CSR_Matrix_T)     :: Miter, LU_Miter  ! compressed row
  TYPE(SpRowColD_T)      :: temp_LU_Dec           ! sparse-LU matrix format
  REAL(dp), ALLOCATABLE  :: rhs(:)

  INTEGER,  ALLOCATABLE  :: LU_Perm(:)

  ! permutation vector/ pivot order for LU-decomp
  INTEGER, ALLOCATABLE :: InvPermu(:)
  INTEGER, ALLOCATABLE :: PivOrder(:)

  ! Timer variables
  REAL(dp) :: StartTimer, TimeSymbolic, TimeNumeric, TimeSolve

  INTEGER :: i

  CHARACTER(400) :: io_msg
  INTEGER        :: io_stat




  !----------------------------------------------------------------
  ! --- Read run control parameters (which runfile)
  CALL getarg( 1 , FileName0 )             
  IF ( FileName0 == '' ) THEN
    WRITE(*,777,ADVANCE='NO') 'Input RUNFilename: '; READ(*,*)   FileName0
  END IF
  FileName0 = TRIM(ADJUSTL(FileName0))

  
  Miter = ReadSparseMatrix(FileName0)

  CALL CPU_TIME(StartTimer)
  !----------------------------------------------------------------
  ! --- convert to Sparse RowColumn format for factorisation
  temp_LU_Dec = CSR_to_SpRowColD(Miter) 

  !----------------------------------------------------------------
  ! --- symbolic factorization

  WRITE(*,*); WRITE(*,777) 'Ordering options, type:'; WRITE(*,*)
  WRITE(*,777) '      marko     ==> Markowitz ordering heuristic (minimum-degree)'
  WRITE(*,777) '      given     ==> pivot order is provided by an integer array'
  WRITE(*,777) '      inorder   ==> pivot order is not altered'
  WRITE(*,*)
  WRITE(*,777,ADVANCE='NO') 'Choose ordering :: '
  READ(*,*,IOSTAT=io_stat,IOMSG=io_msg) input

  IF ( TRIM(input) == 'marko' ) THEN

    CALL SymbLU_SpRowColD_M( temp_LU_Dec )
  
  ELSE IF ( TRIM(input) == 'given' ) THEN
    
    WRITE(*,777,ADVANCE='NO') 'Select File :: '
    READ(*,*,IOSTAT=io_stat,IOMSG=io_msg) input

    PivOrder = InputPivotOrder(input,temp_LU_Dec%m)
    CALL SymbLU_SpRowColD( temp_LU_Dec , PivOrder )
  
  ELSE IF ( TRIM(input) == 'inorder' ) THEN

    PivOrder = [(i , i = 1 , temp_LU_Dec%m )]
    CALL SymbLU_SpRowColD( temp_LU_Dec , PivOrder )

  ELSE

    WRITE(*,777) 'Wrong input --> STOP!'
    STOP

  END IF

  !----------------------------------------------------------------
  ! --- converting back to csr format
  LU_Miter = RowColD_to_CSR( temp_LU_Dec )

  !----------------------------------------------------------------
  ! --- Get the permutation vector LU_Perm and map values of Miter to the permuted LU matrix
  CALL Get_LU_Permutaion( LU_Perm , LU_Miter , Miter )
  !----------------------------------------------------------------
  CALL CPU_TIME(TimeSymbolic)
  
  WRITE(*,*)
  WRITE(*,'(10X,A,1X,F12.8,A)') 'Time symbolic factorization = ', TimeSymbolic - StartTimer, ' [sec]'


  ALLOCATE(rhs(LU_Miter%m))

  CALL RANDOM_NUMBER( LU_Miter%Val )
  CALL RANDOM_NUMBER( rhs )

  CALL CPU_TIME(StartTimer) 
  CALL SparseLU( LU_Miter )
  CALL CPU_TIME(TimeNumeric)
  WRITE(*,'(10X,A,1X,F12.8,A)') 'Time numeric factorization = ', TimeNumeric-StartTimer, ' [sec]'

  CALL CPU_TIME(StartTimer) 
  CALL SolveSparse( LU_Miter , rhs )
  CALL CPU_TIME(TimeSolve)
  WRITE(*,'(10X,A,1X,F12.8,A)') 'Time solving triangula systems = ', TimeSolve-StartTimer, ' [sec]'
  

  !================================================================
  !==  FORMAT Statements
  !================================================================
  !
  777  FORMAT(10X,A)
END PROGRAM Main_SparseGaussElimmination
