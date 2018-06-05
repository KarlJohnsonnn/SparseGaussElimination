!================================================================================!
!                                                                                !
!                         This module includes a collection                      !
!                    of sparse matrix calculations for chemical                  !
!                       reaction systems, main format CSR                        !
!                                                                                !
!================================================================================!
!                                                                                
MODULE Sparse_Mod
  ! Contains:
  !   - SYMBOLIC MATRIX * MATRIC
  !   - SYMBOLIC MATRIX + MATRIX
  !   - TRANSPOSE MATRIX
  !   - SPARSE MATRIX * MATRIC
  !   - SPARSE MATRIX +- MATRIC
  !   - SORT ALGORITHM FOR SYMBOLIC MATRIX*MATRIX CALC
  !   - SPARSE JACOBIAN MATRIX CALC
  !   - SPARSE MITER MATRIX CALC
  !   - CONVERT COMPRESSED ROW FORMAT TO ROW-INDEX, COL-INDEX
  !   - PRINT SPARSE MATRIX (compressed Row format)
  !   - PRINT SPARSE MATRIX (ROW-INDEX COLUMN-INDEX)
  !   - SPARSE MATRIX*VECTOR+VECTOR (SPARSE-MATRIX, DENSE VECTORS)
  ! 
  USE Kind_Mod
  USE mo_unirnk
  !
  IMPLICIT NONE 
  ! 
  INTEGER, PARAMETER, PRIVATE :: inilen=400
  !
  !
  TYPE CSR_Matrix_T            !Compressed Rowindex, standart columnindex
    INTEGER               :: m=0, n=0, nnz=0
    INTEGER, ALLOCATABLE  :: RowPtr(:)
    INTEGER, ALLOCATABLE  :: ColInd(:)
    INTEGER, ALLOCATABLE  :: DiagPtr(:)         ! position of diag entries
    INTEGER, ALLOCATABLE  :: DiagPtr_P(:)       ! permuted pos of diag entries
    INTEGER, ALLOCATABLE  :: DiagPtr_R(:)       ! permuted pos of rate entries
    INTEGER, ALLOCATABLE  :: DiagPtr_C(:)       ! permuted pos of conenctations
    INTEGER, ALLOCATABLE  :: RowVectorPtr(:)    ! for Teq matrix ( -U^T*D_c )
    INTEGER, ALLOCATABLE  :: ColVectorPtr(:)    ! for Teq matrix ( ~K )
    INTEGER               :: XPtr=-42           ! for Teq matrix ( X )
    INTEGER, ALLOCATABLE  :: Permu(:)           ! permutation vector (markowitz)
    INTEGER, ALLOCATABLE  :: InvPer(:)          ! invers permutation (inv markowitz)
    INTEGER, ALLOCATABLE  :: LUperm(:)          ! 
    REAL(dp), ALLOCATABLE :: Val(:)
  END TYPE CSR_Matrix_T
  !
  TYPE SpRowIndColInd_T     !standart Rowindex, standart columnindex
    INTEGER :: m=0,n=0,nnz=0
    INTEGER, ALLOCATABLE :: RowInd(:)
    INTEGER, ALLOCATABLE :: ColInd(:)
    REAL(dp), ALLOCATABLE :: Val(:)
  END TYPE SpRowIndColInd_T
  !
  TYPE SpRowColD_T
    INTEGER :: m=0,n=0
    INTEGER, POINTER :: RowPtr(:,:)=>NULL()
    INTEGER, POINTER :: ColInd(:)=>NULL()
    INTEGER, POINTER :: Permu(:)=>NULL()
    INTEGER, POINTER :: InvPer(:)=>NULL()
    INTEGER :: ep=1
    INTEGER :: last=0
    INTEGER :: len=0
    INTEGER :: nnz=0
  END TYPE SpRowColD_T
  
 
  TYPE(SpRowIndColInd_T) :: MiterFact
  

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE DAX_sparse
    !MODULE PROCEDURE DAX_sparse_n3
  END INTERFACE

  !
  CONTAINS
  !
  !
  !
  FUNCTION New_CSR(m,n,nnz,ri,ci,val) RESULT(newA)
    INTEGER :: m, n
    INTEGER, OPTIONAL :: nnz
    INTEGER, OPTIONAL :: ri(:), ci(:)
    REAL(dp), OPTIONAL :: val(:)
    TYPE(CSR_Matrix_T) :: newA
    !
    INTEGER :: i, j, sameRow, cCnt

    newA%m = m
    newA%n = n
    !
    ALLOCATE(newA%RowPtr(m+1))
    newA%RowPtr = 0
    newA%RowPtr(1) = 1
    !
    IF (PRESENT(nnz)) THEN
      ALLOCATE(newA%ColInd(nnz))
      newA%ColInd = -1
      ALLOCATE(newA%Val(nnz))
      newA%Val = 0.0_dp
      newA%nnz = nnz
    END IF
    
    ! if row and column indices are given
    IF (PRESENT(ri).AND.PRESENT(ci)) THEN
      IF (SIZE(ri) /= SIZE(ci)) STOP ' SIZE(ri) /= SIZE(ci) '
      cCnt = 0
      DO i = 1,m
        sameRow = COUNT(ri==i)
        newA%RowPtr(i+1) = newA%RowPtr(i) + sameRow
      END DO
      newA%ColInd = ci
      IF (PRESENT(val)) newA%val = val
    END IF
  END FUNCTION New_CSR
  !
  !
  SUBROUTINE Free_Matrix_CSR(A)
    TYPE(CSR_Matrix_T) :: A
    !
    IF (ALLOCATED(A%RowPtr))    DEALLOCATE(A%RowPtr)
    IF (ALLOCATED(A%ColInd))    DEALLOCATE(A%ColInd)
    IF (ALLOCATED(A%DiagPtr))   DEALLOCATE(A%DiagPtr)
    IF (ALLOCATED(A%DiagPtr_R)) DEALLOCATE(A%DiagPtr_R)
    IF (ALLOCATED(A%DiagPtr_C)) DEALLOCATE(A%DiagPtr_C)
    IF (ALLOCATED(A%Permu))     DEALLOCATE(A%Permu)
    IF (ALLOCATED(A%InvPer))    DEALLOCATE(A%InvPer)  
    IF (ALLOCATED(A%Val))       DEALLOCATE(A%Val)
  END SUBROUTINE Free_Matrix_CSR
  !
  !
  SUBROUTINE Free_SpRowColD(A)
    TYPE(SpRowColD_T) :: A
    !
    !WRITE(*,*) ' vor free :: '
    !WRITE(*,*) ' m = ', A%m
    !WRITE(*,*) ' n = ', A%n
    !WRITE(*,*) ' ep = ', A%ep
    !WRITE(*,*) ' last = ', A%last
    !WRITE(*,*) ' len = ', A%len
    !WRITE(*,*) ' nnz = ', A%nnz
    !WRITE(*,*) ' size RP = ', SIZE(A%RowPtr)
    !WRITE(*,*) ' size cI = ', SIZE(A%ColInd)
    !WRITE(*,*) ' size P  = ', SIZE(A%Permu)
    !WRITE(*,*) ' size iP = ', SIZE(A%InvPer)

    A%m=0
    A%n=0
    A%ep=1
    A%last=0
    A%len=0
    A%nnz=0
    IF (ASSOCIATED(A%RowPtr)) NULLIFY(A%RowPtr)
    IF (ASSOCIATED(A%ColInd)) NULLIFY(A%ColInd)
    IF (ASSOCIATED(A%Permu )) NULLIFY(A%Permu)
    IF (ASSOCIATED(A%InvPer)) NULLIFY(A%InvPer)
  END SUBROUTINE Free_SpRowColD
  !
  !
  SUBROUTINE Free_SpRowIndColInd(A)
    TYPE(SpRowIndColInd_T) :: A
    !
    IF (ALLOCATED(A%RowInd)) DEALLOCATE(A%RowInd)
    IF (ALLOCATED(A%ColInd)) DEALLOCATE(A%ColInd)
    IF (ALLOCATED(A%Val))    DEALLOCATE(A%Val)
  END SUBROUTINE Free_SpRowIndColInd   
  !
  !
  FUNCTION SparseID(dim) RESULT(Mat)
    TYPE(CSR_Matrix_T) :: Mat
    INTEGER :: dim
    INTEGER :: i
    !
    Mat = New_CSR(dim,dim,dim)
    DO i=1,dim
      Mat%RowPtr(i+1)=Mat%RowPtr(i)+1
      Mat%ColInd(i)=i
    END DO
    Mat%Val = 1.0_dp
  END FUNCTION SparseID
  !

  FUNCTION CSR_to_FULL(CSR) RESULT(Full)
    TYPE(CSR_Matrix_T) :: CSR
    REAL(dp) :: Full(CSR%m,CSR%n)

    INTEGER :: i,j,jj 
    
    Full = 0.0_dp
    DO i=1,CSR%m
      DO jj=CSR%RowPtr(i),CSR%RowPtr(i+1)-1
        j = CSR%ColInd(jj)
        Full(i,j) = CSR%Val(jj)
      END DO
    END DO
  END FUNCTION CSR_to_FULL



  SUBROUTINE CompressIntegerArray(Array)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: Array(:)
    INTEGER, ALLOCATABLE :: tmpArray(:)
    
    INTEGER :: i, N, cnt, M
    
    N = COUNT(Array/=0)
    ALLOCATE(tmpArray(N))

    cnt = 0
    DO i=1,SIZE(Array)
      IF (Array(i)/=0) THEN
        cnt = cnt + 1
        tmpArray(cnt) = Array(i)
      END IF
    END DO
    Array = [tmpArray]
  END SUBROUTINE CompressIntegerArray



  SUBROUTINE CompressDoubleArray(Array)
    REAL(dp), ALLOCATABLE, INTENT(INOUT) :: Array(:)
    REAL(dp), ALLOCATABLE :: tmpArray(:)
    
    INTEGER :: i, N, cnt
    REAL(dp), PARAMETER :: big = -99999999999999.d0
    
    N = COUNT( Array /= big )
    ALLOCATE(tmpArray(N))

    cnt = 0
    DO i=1,SIZE(Array)
      IF ( Array(i) /= big ) THEN
        cnt = cnt + 1
        tmpArray(cnt) = Array(i)
      END IF
    END DO
    Array = [tmpArray]
  END SUBROUTINE CompressDoubleArray


  FUNCTION FULL_to_CSR(Full) RESULT(CSR)
    REAL(dp) :: Full(:,:)
    TYPE(CSR_Matrix_T) :: CSR

    INTEGER :: i,j,jj 
    INTEGER :: m,n,nnz 

    m = SIZE(Full,1)
    n = SIZE(Full,2)
    nnz = COUNT(Full /= 0.0_dp)

    CSR = New_CSR(m,n,nnz)
    
    jj = 0
    DO i = 1 , m
      CSR%RowPtr(i+1) = CSR%RowPtr(i) + COUNT(Full(i,:)/=0.0_dp)
      DO j = 1 , n
        IF ( Full(i,j) /= 0.0_dp) THEN
          jj = jj + 1
          CSR%ColInd(jj) = j
          CSR%Val(jj)    = Full(i,j)
        END IF
      END DO
    END DO
    CSR%nnz = jj
  END FUNCTION FULL_to_CSR


  FUNCTION Copy_SpRowIndColInd(orig) RESULT(copy)
    TYPE(SpRowIndColInd_T) :: orig
    TYPE(SpRowIndColInd_T) :: copy

    copy%m = orig%m
    copy%n = orig%n
    copy%nnz = orig%nnz
    copy%RowInd = orig%RowInd
    copy%ColInd = orig%ColInd
    copy%Val    = orig%Val
  END FUNCTION Copy_SpRowIndColInd


  FUNCTION Copy_SpRowColD(orig) RESULT(copy)
    TYPE(SpRowColD_T) :: orig
    TYPE(SpRowColD_T) :: copy

    !ALLOCATE(copy%RowInd()
    copy%m = orig%m
    copy%n = orig%n
    copy%nnz = orig%nnz
    copy%RowPtr = orig%RowPtr
    copy%ColInd = orig%ColInd
    copy%Permu  = orig%Permu
    copy%InvPer = orig%InvPer

    !INTEGER :: m=0,n=0
    !INTEGER, POINTER :: RowPtr(:,:)=>NULL()
    !INTEGER, POINTER :: ColInd(:)=>NULL()
    !INTEGER, POINTER :: Permu(:)=>NULL()
    !INTEGER, POINTER :: InvPer(:)=>NULL()
    !INTEGER :: ep=1
    !INTEGER :: last=0
    !INTEGER :: len=0
    !INTEGER :: nnz=0
  END FUNCTION Copy_SpRowColD


  FUNCTION Copy_CSR(orig) RESULT(copy)
    TYPE(CSR_Matrix_T) :: orig
    TYPE(CSR_Matrix_T) :: copy

    copy = New_CSR(orig%m,orig%n,orig%nnz)
    copy%RowPtr = orig%RowPtr
    copy%ColInd = orig%ColInd
    copy%Val    = orig%Val
    Copy%nnz    = orig%nnz
    IF (orig%XPtr/=-42) copy%XPtr = orig%XPtr
    IF (ALLOCATED(orig%DiagPtr))   copy%DiagPtr   = orig%DiagPtr
    IF (ALLOCATED(orig%DiagPtr_P)) copy%DiagPtr_P = orig%DiagPtr_P
    IF (ALLOCATED(orig%DiagPtr_R)) copy%DiagPtr_R = orig%DiagPtr_R
    IF (ALLOCATED(orig%DiagPtr_C)) copy%DiagPtr_C = orig%DiagPtr_C
    IF (ALLOCATED(orig%RowVectorPtr)) copy%RowVectorPtr = orig%RowVectorPtr
    IF (ALLOCATED(orig%ColVectorPtr)) copy%ColVectorPtr = orig%ColVectorPtr
    IF (ALLOCATED(orig%Permu))  copy%Permu  = orig%Permu
    IF (ALLOCATED(orig%InvPer)) copy%InvPer = orig%InvPer
    IF (ALLOCATED(orig%LUperm)) copy%LUperm = orig%LUperm
  END FUNCTION Copy_CSR

 
  FUNCTION RowColD_to_CSR(SpRow) RESULT(CSR)
    TYPE(CSR_Matrix_T) :: CSR
    TYPE(SpRowColD_T)  :: SpRow
    !
    INTEGER :: i,j,jj,nzrA
    ! 

    CSR = New_CSR(SpRow%m,SpRow%n,SpRow%nnz)
    
    ALLOCATE(CSR%DiagPtr(SpRow%m),CSR%DiagPtr_P(SpRow%m))
  
    CSR%DiagPtr   = -11
    CSR%DiagPtr_P = -11
    
    NzrA=0
    DO i=1,CSR%m
      CSR%RowPtr(i+1)=CSR%RowPtr(i)
      DO jj=SpRow%RowPtr(1,i),SpRow%RowPtr(2,i)
        j=SpRow%ColInd(jj) 
        NzrA=NzrA+1
        CSR%ColInd(NzrA)=j
        IF (i==j) THEN
          CSR%DiagPtr(i)=CSR%RowPtr(i+1)
        END IF
        CSR%RowPtr(i+1)=CSR%RowPtr(i+1)+1
      END DO
    END DO
    
    IF (ASSOCIATED(SpRow%Permu)) THEN
      IF (.NOT.ALLOCATED(CSR%Permu)) THEN
        ALLOCATE(CSR%Permu(CSR%n))
        CSR%Permu(:)=-16
      END IF
      CSR%Permu(:)=SpRow%Permu(:)
    END IF
    IF (ASSOCIATED(SpRow%InvPer)) THEN
      IF (.NOT.ALLOCATED(CSR%InvPer)) THEN
        ALLOCATE(CSR%InvPer(CSR%n))
        CSR%InvPer(:)=-16
      END IF
      CSR%InvPer(:)=SpRow%InvPer(:)
    END IF

  END FUNCTION RowColD_to_CSR



  FUNCTION SpRowIndColInd_to_CSR(SpRiCi) RESULT(CSR)
    ! IN:
    TYPE(SpRowIndColInd_T) :: SpRiCi
    ! OUT:
    TYPE(CSR_Matrix_T)     :: CSR
    ! TEMP:
    INTEGER :: i, ii, jj, n, m, nnz

    m = SpRiCi%m
    n = SpRiCi%n
    nnz = SpRiCi%nnz

    CSR = New_CSR(m,n,nnz)
    
    ii = 0
    DO i = 1 , m
      CSR%RowPtr(i+1) = CSR%RowPtr(i) + COUNT(SpRiCi%RowInd==i)
      DO jj = CSR%RowPtr(i) , CSR%RowPtr(i+1)-1
        ii = ii + 1
        CSR%ColInd(jj) = SpRiCi%ColInd(ii)
        CSR%Val(jj)    = SpRiCi%Val(ii)
      END DO
    END DO
    CSR%nnz = ii
  END FUNCTION SpRowIndColInd_to_CSR
  


  !
  FUNCTION CSR_to_SpRowColD(CSR) RESULT(SpRowCol)
    TYPE(SpRowColD_T) :: SpRowCol
    TYPE(CSR_Matrix_T) :: CSR
    !
    INTEGER :: nzrCSR
    INTEGER :: AddLen
    INTEGER :: Start,End
    INTEGER :: i,j,jj
    !
    AddLen=10
    SpRowCol%m=CSR%m
    SpRowCol%n=CSR%n
    ALLOCATE(SpRowCol%RowPtr(2,SpRowCol%n))
    nzrCSR=SIZE(CSR%ColInd)
    SpRowCol%len=10*nzrCSR+AddLen*SpRowCol%n
    ALLOCATE(SpRowCol%ColInd(SpRowCol%len))
    SpRowCol%ColInd=0
    ALLOCATE(SpRowCol%Permu(SpRowCol%n))
    SpRowCol%Permu=0
    ALLOCATE(SpRowCol%InvPer(SpRowCol%n))
    SpRowCol%InvPer=0
    !
    Start=1
    DO i=1,CSR%n
      SpRowCol%RowPtr(1,i)=Start
      End=Start+(CSR%RowPtr(i+1)-CSR%RowPtr(i)-1)
      SpRowCol%RowPtr(2,i)=End
      DO jj=CSR%RowPtr(i),CSR%RowPtr(i+1)-1
        j=CSR%ColInd(jj)
        SpRowCol%ColInd(Start)=j
        Start=Start+1
      END DO  
      Start=Start+AddLen-1
    END DO  
    SpRowCol%ep=SpRowCol%RowPtr(2,SpRowCol%n)+1
    SpRowCol%last=SpRowCol%n  
  END FUNCTION CSR_to_SpRowColD
  


  !
  SUBROUTINE Sort_SpRowColD(A)
    TYPE(SpRowColD_T) :: A
    INTEGER :: i,jj
    DO i=1,A%m
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        A%ColInd(jj)=A%Permu(A%ColInd(jj))
      END DO
      CALL SortVec(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
    END DO
  END SUBROUTINE Sort_SpRowColD
  


  !
  SUBROUTINE SymbLU_SpRowColD(A,Permu)
    TYPE(SpRowColD_T) :: A
    INTEGER :: Permu(:)
    !
    INTEGER :: RowPiv(A%n)
    INTEGER :: i,j,l,jj,ip,iPiv
    LOGICAL :: ins
    !
    DO i=1,A%n
      A%InvPer(i)=i
      A%Permu(i)=i
    END DO
    !
    DO i=1,A%n
      ip=A%Permu(Permu(i))
      CALL Swap(A%InvPer(i),A%InvPer(ip))
      CALL Swap(A%RowPtr(1,i),A%RowPtr(1,ip))
      CALL Swap(A%RowPtr(2,i),A%RowPtr(2,ip))
      A%Permu(A%InvPer(i))=i
      A%Permu(A%InvPer(ip))=ip
      IF (A%last==i) THEN
        A%last=ip
      ELSE IF (A%last==ip) THEN
        A%last=i
      END IF
      !
      !   Update
      iPiv=0
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        IF (A%Permu(A%ColInd(jj))>i) THEN
          iPiv=iPiv+1 
          RowPiv(iPiv)=A%ColInd(jj)
        END IF
      END DO
      IF (iPiv>0) THEN
        DO j=i+1,A%n
          DO jj=A%RowPtr(1,j),A%RowPtr(2,j)
            IF (A%Permu(A%ColInd(jj))==i) THEN
              DO l=1,iPiv
                CALL Insert_SpRowColD(A,j,RowPiv(l),ins)
              END DO
              EXIT
            END IF
          END DO
        END DO
      END IF
    END DO
    DO i=1,A%n
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        A%ColInd(jj)=A%Permu(A%ColInd(jj))
      END DO
      CALL SortVec(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
      A%nnz=A%nnz+SIZE(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
    END DO
  END SUBROUTINE SymbLU_SpRowColD
  


  !
  SUBROUTINE SymbLU_SpRowColD_M(A)
    TYPE(SpRowColD_T) :: A
    !
    INTEGER :: r(A%n),c(A%n),RowPiv(A%n)
    INTEGER :: i,j,l,jj,ip,ip1(1),iPiv
    REAL(dp) :: md 
    INTEGER ::  rc
    INTEGER, ALLOCATABLE :: rc_idx(:,:)
    LOGICAL :: ins
    INTEGER :: indexP(1,1)
  
    c = 0
    DO i = 1 , A%n
      A%InvPer(i) = i
      A%Permu(i)  = i
      ! Compute initial Markowitz count
      r(i) = A%RowPtr(2,i) - A%RowPtr(1,i) + 1
      DO jj = A%RowPtr(1,i) , A%RowPtr(2,i)
        c(A%ColInd(jj)) = c(A%ColInd(jj)) + 1
      END DO
    END DO

    !ALLOCATE(rc_idx(2,3))
    !rc_idx = 0

    
    MAIN_LOOP: DO i = 1 , A%n
      ip = 0
      md = 1.d99
      
      !ALLOCATE(rc_idx(0))
      DO j = i , A%n
        rc = (r(j)-1) * (c(j)-1)
        IF ( rc <= md ) THEN
          !rc_idx = [ j , rc_idx]
          md = rc
          ip = j
        END IF
      END DO
      !ip = rc_idx(1)
      !DEALLOCATE(rc_idx)

      !WRITE(*,*) 'ip loop/vect = ', ip, rc_idx(:)
      !READ(*,*)
      
      !! wird nie erreicht?
      !IF (ip==0) THEN
      !  ip1(:) = MINLOC((r(i:A%n)-1)*(c(i:A%n)-1))+(i-1)
      !  ip = ip1(1)
      !  MarkowitzCounts(i,1) = ip1(1)
      !  MarkowitzCounts(i,2) = SIZE(ip1)
      !END IF
      
      CALL Swap( r(i) , r(ip) )
      CALL Swap( c(i) , c(ip) )
      CALL Swap( A%InvPer(i)   , A%InvPer(ip) )
      CALL Swap( A%RowPtr(1,i) , A%RowPtr(1,ip) )
      CALL Swap( A%RowPtr(2,i) , A%RowPtr(2,ip) )
      
      A%Permu(A%InvPer(i))  = i
      A%Permu(A%InvPer(ip)) = ip

      IF ( A%last == i ) THEN
        A%last = ip
      ELSE IF ( A%last == ip ) THEN
        A%last = i
      END IF
      !
      ! Update
      iPiv = 0
      DO jj = A%RowPtr(1,i) , A%RowPtr(2,i)
        IF ( A%Permu(A%ColInd(jj)) > i ) THEN
          iPiv = iPiv + 1 
          RowPiv(iPiv) = A%ColInd(jj)
          c(A%Permu(A%ColInd(jj))) = c(A%Permu(A%ColInd(jj))) - 1
        END IF
      END DO
      IF ( iPiv > 0 ) THEN
        DO j = i+1 , A%n
          DO jj = A%RowPtr(1,j) , A%RowPtr(2,j)
            IF ( A%Permu(A%ColInd(jj)) == i ) THEN
              r(j) = r(j) - 1 
              c(i) = c(i) - 1
              DO l = 1 , iPiv
                CALL Insert_SpRowColD( A, j, RowPiv(l), ins)
                IF (ins) THEN
                  c(A%Permu(RowPiv(l))) = c(A%Permu(RowPiv(l))) + 1
                  r(j) = r(j) + 1
                END IF
              END DO
              EXIT
            END IF
          END DO
          IF ( c(i) == 1 ) EXIT
        END DO
      END IF
    END DO MAIN_LOOP

    DO i = 1 , A%n
      DO jj = A%RowPtr(1,i) , A%RowPtr(2,i)
        A%ColInd(jj) = A%Permu(A%ColInd(jj))
      END DO
      CALL SortVec( A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)) )
      A%nnz = A%nnz + SIZE(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
    END DO
  END SUBROUTINE SymbLU_SpRowColD_M
  !
  !
  SUBROUTINE Swap(i,j)
    INTEGER :: i, j, iTemp
    iTemp=i;    i=j;    j=iTemp
  END SUBROUTINE Swap
  !
  !
  SUBROUTINE SortVec(vec)
    INTEGER :: vec(:)
    !
    INTEGER :: i,itemp,j,n
    n=SIZE(Vec)
    DO i=1,n
      DO j=1,n-i
        IF (vec(j)>vec(j+1)) THEN
          itemp=vec(j)
          vec(j)=vec(j+1)
          vec(j+1)=itemp
        END IF
      END DO
    END DO
  END SUBROUTINE SortVec
  !
  !
  SUBROUTINE SortVecDesc(vec)
    INTEGER :: vec(:)
    !
    INTEGER :: i,itemp,j,n
   
    n = SIZE(Vec)
   
    
    DO i = 1 , n
      DO j = 1 , n-i
        IF ( vec(j) < vec(j+1) ) THEN
          itemp    = vec(j)
          vec(j)   = vec(j+1)
          vec(j+1) = itemp
        END IF
      END DO
    END DO
  END SUBROUTINE SortVecDesc
  !

  SUBROUTINE SortVecDesc2(vec,q)
    INTEGER, INTENT(INOUT) :: vec(:)
    INTEGER, ALLOCATABLE, OPTIONAL :: q(:)
    !
    INTEGER :: i,n,iMax(1)
    INTEGER, ALLOCATABLE :: tmpVec(:)
   
    n = SIZE(Vec)
    ALLOCATE(tmpVec(n));  tmpVec = 0
   
    IF (PRESENT(q).AND..NOT.ALLOCATED(q)) ALLOCATE(q(n))
    
    DO i = 1 , n
      iMax = MAXLOC(vec)
      tmpVec(i) = vec(iMax(1))
      vec(iMax(1)) = -1
      IF (PRESENT(q)) q(i) = iMax(1)
    END DO
    vec = tmpVec
  END SUBROUTINE SortVecDesc2


  !
  SUBROUTINE Insert_SpRowColD(A,iA,jA,ins)
    TYPE(SpRowColD_T) :: A
    INTEGER :: iA,jA
    LOGICAL :: ins
    !
    INTEGER :: itemp,j,l
    INTEGER, ALLOCATABLE :: iWork(:)
    !
    ! Test ob Element (ia,ja) bereits enthalten
    IF (iA==0.OR.jA==0) THEN
      WRITE(*,*) 'iA',iA
      WRITE(*,*) 'jA',jA
      STOP 'STOP'
    END IF  
    ins=.TRUE.
    DO j=A%RowPtr(1,iA),A%RowPtr(2,iA)
      IF (jA==A%ColInd(j)) THEN
        ins=.FALSE.
      END IF
    END DO
    !
    ! Test auf freien Speicherplatz in der ia-ten
    ! Zeile von a
    ! 
    IF (ins) THEN
      IF (A%ColInd(A%RowPtr(2,iA)+1)/=0) THEN
        ! ja-te Zeile von a wird nach hinten
        itemp=A%ep
        DO l=A%RowPtr(1,iA),A%RowPtr(2,iA)
          A%ColInd(A%ep)=A%ColInd(l)
          A%ColInd(l)=0
          A%ep=A%ep+1
        END DO
        A%RowPtr(2,iA)=A%ep-1
        A%RowPtr(1,iA)=itemp
        ! A%ep=A%ep+1
        A%last=iA
      ENDIF
      A%RowPtr(2,iA)=A%RowPtr(2,iA)+1
      A%ColInd(A%RowPtr(2,iA))=jA
      IF (iA==A%last) A%ep=A%ep+1
    END IF
    !
    IF (A%ep>=A%len-A%m) THEN
      CALL gcmat_SpRowColD(A)
    END IF
    !
    IF (A%ep>=A%len-A%m) THEN
      !   Speicherplatz von A nicht ausreichend
      ALLOCATE(iWork(A%ep))
      iWork(1:A%ep)=A%ColInd(1:A%ep)
      DEALLOCATE(A%ColInd)
      A%len=2*A%len
      ALLOCATE(A%ColInd(A%len))
      A%ColInd(1:A%ep)=iWork(1:A%ep)
      DEALLOCATE(iWork)  
    END IF
  END SUBROUTINE Insert_SpRowColD
  !
  ! 
  SUBROUTINE gcmat_SpRowColD(A)
   !   Externe Variable
   TYPE (SpRowColD_T) :: A
   ! 
   !   gcmat komprimiert eine zeilenorientierte, dynamische
   !   Speicherstruktur einer schwachbesetzten Matrix a
   !   der Ordnung n. Die Spaltenindizes der Nichtnull-
   !   elemente der i-ten Zeile von a sind in A%ColInd(A%RowPtr(i,1)),
   !   A%ColInd(A%RowPtr(1,i)+1)...,A%ColInd(A%RowPtr(2,i)) ent-
   !   halten.
   ! 
   !    Beschreibung der Parameter
   ! 
   ! 
   !    A%n      (i/o) integer
   !                 Dimension of matrix a
   !
   !    A%RowPtr (i/o) integer(2,n)
   !                 Pointerfeld zur Beschreibung von a. A%RowPtr muss
   !                 durch das rufende Programm belegt werden.
   ! 
   !    A%ColInd (i/o) integer(len)
   !   Externe Variable
   !   TYPE (SpRowColD_T) :: A
   ! 
   !   gcmat komprimiert eine zeilenorientierte, dynamische
   !   Speicherstruktur einer schwachbesetzten Matrix a
   !   der Ordnung n. Die Spaltenindizes der Nichtnull-
   !   elemente der i-ten Zeile von a sind in A%ColInd(A%RowPtr(i,1)),
   !   A%ColInd(A%RowPtr(1,i)+1)...,A%ColInd(A%RowPtr(2,i)) ent-
   !   halten.
   ! 
   !    Beschreibung der Parameter
   ! 
   ! 
   !    A%n      (i/o) integer
   !                 Dimension of matrix a
   !
   !    A%RowPtr (i/o) integer(2,n)
   !                 Pointerfeld zur Beschreibung von a. A%RowPtr muss
   !                 durch das rufende Programm belegt werden.
   ! 
   !    A%ColInd (i/o) integer(len)
   !                 Feld zur dynamischen Verwaltung der Spaltenindizes
   !                 der Nichtnullelemente von a. 
   ! 
   !    A%ep     (i/o) integer
   !                 Pointer der auf das erste freie Feld in a verweist,
   !                 d.h A%ColInd(ep),...,A%ColInd(len) sind frei verfuegbar.
   ! 
   !    A%len    (i)   integer
   !                 Gesamtlaenge des Feldes A%ColInd. 
   ! 
   !  Interne Variable
   !
   INTEGER i,iz,j,l,pointr,rowlen,ep,len,m
   !
   m=A%m
   ep=A%ep
   len=A%len
   pointr=1
   i=1
   !
   DO 
     IF (i>=ep) EXIT
       IF (A%ColInd(i).ne.0) THEN
         !
         ! Ermittlung der aktuellen Zeile sowie deren Laenge
         DO l=1,m
           IF (A%RowPtr(1,l).le.i.and.i.le.A%RowPtr(2,l)) THEN
             iz=l
           END IF
         END DO
         rowlen=A%RowPtr(2,iz)-A%RowPtr(1,iz)
         ! 
         ! Setzen der neuen Anfangs- und Endadresse der
         ! aktuellen Zeile
         A%RowPtr(1,iz)=pointr
         A%RowPtr(2,iz)=pointr+rowlen
         DO j=pointr,pointr+rowlen
           A%ColInd(j)=A%ColInd(i)
           i=i+1
         END DO
         i=i-1
         pointr=A%RowPtr(2,iz)+1
       ENDIF
      i=i+1
    END DO
    !
    !   Belegung des freien Teils von A%ColInd mit 0
    !
    ep=pointr
    DO i=1,m
      IF (A%RowPtr(1,i).gt.ep) THEN
        A%RowPtr(1,i)=ep
        A%RowPtr(2,i)=A%RowPtr(1,i)-1
        ep=ep+inilen
      END IF
    END DO
    !        
    DO i=pointr,len
      A%ColInd(i)=0
    END DO
    A%ep=ep
  END SUBROUTINE gcmat_SpRowColD
  !
  !
  SUBROUTINE SpDeallocate_SpRowColD(A)
    TYPE(SpRowColD_T) :: A
    !
    A%m=0
    A%n=0
    IF (ASSOCIATED(A%RowPtr)) DEALLOCATE(A%RowPtr)
    IF (ASSOCIATED(A%ColInd)) DEALLOCATE(A%ColInd)
    IF (ASSOCIATED(A%Permu))  DEALLOCATE(A%Permu)
    IF (ASSOCIATED(A%InvPer)) DEALLOCATE(A%InvPer)
  END SUBROUTINE SpDeallocate_SpRowColD
  



  ! Transpose Matrix
  SUBROUTINE TransposeSparse(MatAT,MatA)
    TYPE(CSR_Matrix_T), INTENT(OUT) :: MatAT
    TYPE(CSR_Matrix_T), INTENT(IN)  :: MatA
    !
    INTEGER :: i,j                                 ! zählvariabl schleIFen
    INTEGER :: indx
    !
    !
    MatAT = New_CSR(MatA%n,MatA%m)
    !
    !
    DO i=1,MatA%m
      DO j=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        MatAT%RowPtr(MatA%ColInd(j)+1)=MatAT%RowPtr(MatA%ColInd(j)+1)+1
      END DO
    END DO
    !
    DO i=1,MatAT%m
      !print*, 'RowPtr=',MatAT%RowPtr(i)
      MatAT%RowPtr(i+1)=MatAT%RowPtr(i)+MatAT%RowPtr(i+1)
    END DO
      !print*, 'RowPtr=',MatAT%RowPtr(MatAT%m+1)
    !
    ALLOCATE(MatAT%ColInd(SIZE(MatA%ColInd)))
    ALLOCATE(MatAT%Val(SIZE(MatA%Val)))
    MatAT%ColInd=0
    MatAT%Val=0.0_dp
    DO i=1,MatA%m
      DO j=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        indx=MatA%ColInd(j)
        MatAT%ColInd(MatAT%RowPtr(indx))=i
        MatAT%Val(MatAT%RowPtr(indx))=MatA%Val(j)
        MatAT%RowPtr(indx)=MatAT%RowPtr(indx)+1
      END DO
    END DO
    DO i=MatAT%m,1,-1
      MatAT%RowPtr(i+1)=MatAT%RowPtr(i)
    END DO
    MatAT%RowPtr(1)=1
    MatAT%nnz=MatA%RowPtr(MatA%m+1)-1
  END SUBROUTINE TransposeSparse
  


  
  ! SYMBOLIC MATRIX * MATRIX
  SUBROUTINE SymbolicMult(A,B,C)
    ! A*B=C
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    TYPE(CSR_Matrix_T), INTENT(IN) :: B
    TYPE(CSR_Matrix_T), INTENT(OUT) :: C
    !
    INTEGER ::  indx(MAX(A%n,A%m,B%n))
    INTEGER :: i, j, jj, k
    INTEGER :: istart, length, iTemp
    !
      !symbolic matrix multiply c=a*b
    C = New_CSR(A%m,B%n)
    !
    !main loop
    indx=0
    DO i=1,A%m
      iStart=-1
      Length=0
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        DO k=B%RowPtr(j),B%RowPtr(j+1)-1
          IF (indx(B%ColInd(k))==0) THEN
            indx(B%ColInd(k))=istart
            istart=B%ColInd(k)
            length=length+1
          END IF
        END DO
      END DO
      C%RowPtr(i+1)=C%RowPtr(i)+Length
      !
      DO j=C%RowPtr(i),C%RowPtr(i+1)-1
        iTemp=iStart
        istart=indx(istart)
        indx(iTemp)=0
      END DO
      indx(i) = 0
    END DO
    !==========================================
    ALLOCATE(C%ColInd(C%RowPtr(C%m+1)-1))
    C%ColInd=0
    !
    DO i=1,A%m
      iStart=-1
      Length=0
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        DO k=B%RowPtr(j),B%RowPtr(j+1)-1
          IF (indx(B%ColInd(k))==0) THEN
            indx(B%ColInd(k))=istart
            istart=B%ColInd(k)
            length=length+1
          END IF
        END DO
      END DO
      C%RowPtr(i+1)=C%RowPtr(i)+length
      DO j=C%RowPtr(i),C%RowPtr(i+1)-1
       C%ColInd(j)=istart
       istart=indx(istart)
       indx(C%ColInd(j))=0
      END DO
      indx(i) = 0
    END DO
    !
    DO i=1,C%m
      CALL Sort(C%Colind(C%RowPtr(i):C%RowPtr(i+1)-1))
    END DO
    ALLOCATE(C%Val(C%RowPtr(C%m+1)-1))
    C%Val=0.0_dp
    C%nnz=C%RowPtr(C%m+1)-1
  END SUBROUTINE SymbolicMult
  
  


  ! SPARSE MATRIX*MATRIX CALC
  SUBROUTINE SparseMult(A,B,C)
    ! A * B = C
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: A
    TYPE(CSR_Matrix_T), INTENT(IN) :: B
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: C
    !
    INTEGER :: i,j,jj,k,kk
    REAL(dp) :: ajj
    REAL(dp) :: temp(MAX(A%m,A%n,B%n))
    temp=0.0_dp
    !
    DO i=1,A%m
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        ajj=A%Val(jj)
        DO kk=B%RowPtr(j),B%RowPtr(j+1)-1
          k=B%ColInd(kk)
          temp(k)=temp(k)+ajj*B%Val(kk)
        END DO
      END DO
      DO j=C%RowPtr(i),C%RowPtr(i+1)-1
        C%Val(j)=temp(C%ColInd(j))
        temp(C%ColInd(j))=0.
      END DO
    END DO
  END SUBROUTINE SparseMult



 
  ! SYMBOLIC MATRIX + MATRIX  -->  MatC = MatA + MatB
  SUBROUTINE SymbolicAdd(MatC,MatA,MatB)
    TYPE(CSR_Matrix_T), INTENT(IN)  :: MatA
    TYPE(CSR_Matrix_T), INTENT(IN)  :: MatB
    TYPE(CSR_Matrix_T), INTENT(OUT) :: MatC
    !
    INTEGER :: i,ii,j,jj,k,kk
    INTEGER, ALLOCATABLE :: TmpCol(:)
    INTEGER, ALLOCATABLE :: PermVec(:)
    INTEGER, ALLOCATABLE :: Indizes(:)
    INTEGER :: ColLen
    INTEGER :: currentlength
    INTEGER :: sameCnt
    !
    IF (.NOT.((MatA%m==MatB%m).AND.(MatA%n==MatB%n))) THEN
      WRITE(*,*) 'Wrong Matrix Dim'
      WRITE(*,*) 'A: ',MatA%m,MatA%n
      WRITE(*,*) 'B: ',MatB%m,MatB%n
      STOP 'STOP'
    END IF
    !
    MatC = New_CSR(MatA%m,MatA%n)
    !
    DO i=1,MatC%m
      currentlength=(MatA%RowPtr(i+1)-MatA%RowPtr(i))+(MatB%RowPtr(i+1)-MatB%RowPtr(i))
      ALLOCATE(Indizes(currentlength))
      Indizes=0
      sameCnt=0
      DO ii=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        DO j=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
          IF (MatA%ColInd(ii)==MatB%ColInd(j)) THEN
            sameCnt=sameCnt+1
            Indizes(sameCnt)=j
          END IF
        END DO
      END DO
      MatC%RowPtr(i+1)=MatC%RowPtr(i)+currentlength-sameCnt
      DEALLOCATE(Indizes)
    END DO
    !
    ! Allocate ColInd
    ALLOCATE(MatC%ColInd(MatC%RowPtr(MatC%m+1)-1))
    MatC%ColInd=0
    !
    kk=1
    DO i=1,MatC%m
      k=1
      currentlength=MatC%RowPtr(i+1)-MatC%RowPtr(i)
      sameCnt=0
      DO ii=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        DO j=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
          IF (MatA%ColInd(ii)==MatB%ColInd(j)) THEN
            sameCnt=sameCnt+1
          END IF
        END DO
      END DO
      !
      ALLOCATE(TmpCol(currentlength+sameCnt))
      TmpCol=0
      ALLOCATE(PermVec(currentlength+sameCnt))
      PermVec=0
      !
      DO jj=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        TmpCol(k)=MatA%ColInd(jj)
        k=k+1
      END DO
      DO jj=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
        TmpCol(k)=MatB%ColInd(jj)
        k=k+1
      END DO
      !
      CALL unirnk(TmpCol,PermVec,ColLen)
      DO j=1,ColLen
        MatC%ColInd(kk)=TmpCol(PermVec(j))
        kk=kk+1
      END DO
      !
      DEALLOCATE(TmpCol)
      DEALLOCATE(PermVec)
    END DO
    !
    ALLOCATE(MatC%Val(MatC%RowPtr(MatC%m+1)-1))
    MatC%Val=0.0_dp
    MatC%nnz=MatC%RowPtr(MatC%m+1)-1
  END SUBROUTINE SymbolicAdd
  !
  !
  ! SPARSE  MATRIX + MATRIX  -->  MatC = MatA (Sub_Add) MatB 
  SUBROUTINE SparseAdd(MatC,MatA,MatB,Sub)
    TYPE(CSR_Matrix_T), INTENT(IN)    :: MatA, MatB
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: MatC
    CHARACTER, OPTIONAL,  INTENT(IN)    :: Sub
    !
    ! Temp variables
    INTEGER, ALLOCATABLE :: tmpCol(:)
    INTEGER, ALLOCATABLE :: perm(:)
    REAL(dp), ALLOCATABLE :: tmpBVal(:)
    REAL(dp), ALLOCATABLE :: tmpVal(:)
    INTEGER :: lenColIndA, lenColIndB
    INTEGER :: i, ii, j, jj, k ,kk
    !
    ALLOCATE(tmpBVal(MatB%RowPtr(MatB%m+1)-1))
    IF (Sub=='-') THEN
      tmpBVal = -MatB%Val
    ELSE
      tmpBVal = MatB%Val
    END IF
    !
    DO i=1,MatC%m
      DO ii=MatC%RowPtr(i),MatC%RowPtr(i+1)-1
        DO jj=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
          IF (MatA%ColInd(jj)==MatC%ColInd(ii)) THEN
            MatC%Val(ii) = MatC%Val(ii)+MatA%Val(jj)
          END IF
        END DO
        DO kk=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
          IF (MatB%ColInd(kk)==MatC%ColInd(ii)) THEN
            MatC%Val(ii) = MatC%Val(ii)+tmpBVal(kk)
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(tmpBVal)
  END SUBROUTINE SparseAdd
  
  

  ! SORT ALGORITHM FOR SYMBOLIC MATRIX*MATRIX CALC
  SUBROUTINE Sort(v)
    INTEGER :: v(:)
    !
    INTEGER :: i,j,temp
    DO i=1,SIZE(v)
      DO j=i+1,SIZE(v)
        IF (v(i)>v(j)) THEN
          temp=v(i)
          v(i)=v(j)
          v(j)=temp
        END IF
      END DO
    END DO
  END SUBROUTINE Sort
  
  

  ! SORT ALGORITHM FOR MATRIX+MATRIX
  SUBROUTINE Sort2(v,perm)
    INTEGER :: v(:)
    INTEGER :: perm(:)
    !
    INTEGER :: i,j,temp1,temp2
    DO i=1,SIZE(v)
      perm(i)=i
    END DO
    !
    DO i=1,SIZE(v)
      DO j=i+1,SIZE(v)
        IF (v(i)>v(j)) THEN
          ! Feldelemente sortieren
          temp1=v(i)
          v(i)=v(j)
          v(j)=temp1
          ! Permutationsvektor erzeugen
          temp2=perm(i)
          perm(i)=perm(j)
          perm(j)=temp2
        END IF
      END DO
    END DO
  END SUBROUTINE Sort2
  
  
  
  SUBROUTINE CompressList(ColInd,Val)
    INTEGER, ALLOCATABLE :: ColInd(:)
    REAL(dp), ALLOCATABLE, OPTIONAL :: Val(:)
    !
    INTEGER :: i,j,iList,MemberCol
    INTEGER :: TempListCol(SIZE(ColInd))
    REAL(dp) :: MemberVal
    REAL(dp) :: TempListVal(SIZE(ColInd))
    LOGICAL :: Insert
    !
    TempListVal=0.0_dp
    iList=0
    !
    S1:DO i=1,SIZE(ColInd)
      MemberCol=ColInd(i)
      IF (PRESENT(Val)) MemberVal=Val(i)
      Insert=.TRUE.
      S2:DO j=1,iList
        IF (MemberCol==TempListCol(j)) THEN
          Insert=.FALSE.
          IF (PRESENT(Val)) TempListVal(iList)=TempListVal(iList)+MemberVal
          EXIT S2
        END IF
      END DO S2
      IF (Insert) THEN
        iList=iList+1
        TempListCol(iList)=MemberCol
        IF (PRESENT(Val)) TempListVal(iList)=MemberVal
      END IF
    END DO S1
    DEALLOCATE(ColInd)
    IF (PRESENT(Val)) DEALLOCATE(Val)
    ALLOCATE(ColInd(1:iList))
    IF (PRESENT(Val)) ALLOCATE(Val(1:iList))
    ColInd=TempListCol(1:iList)
    IF (PRESENT(Val)) Val=TempListVal(1:iList)
  END SUBROUTINE CompressList




  ! Permutes the values in Miter and writes it to LU structur of Miter
  ! Permutation vector is generated in this routine
  SUBROUTINE Get_LU_Permutaion(Permutation,LU,A)
    ! INOUT
    TYPE(CSR_Matrix_T) :: LU
    ! IN
    TYPE(CSR_Matrix_T) :: A
    INTEGER            :: nnzA
    ! OUT
    INTEGER, ALLOCATABLE :: Permutation(:)
    ! TEMP
    INTEGER :: i,ip,j,jj,jp,jjp,jp1
 
    nnzA=A%RowPtr(A%m+1)-1
   
    IF (.NOT.ALLOCATED(Permutation)) ALLOCATE(Permutation(nnzA))
    Permutation = -14

    LU%Val = 0.0_dp

    DO i = 1 , A%n  
      ip = LU%Permu(i)
      DO jj = A%RowPtr(i) , A%RowPtr(i+1)-1
        j  = A%ColInd(jj)
        jp = LU%Permu(A%ColInd(jj))
        DO jjp = LU%RowPtr(ip) , LU%RowPtr(ip+1)-1
          jp1 = LU%ColInd(jjP)
          IF ( jp1 == jp ) THEN
            LU%Val(jjP) = A%Val(jj)
            Permutation(jj) = jjP
          END IF  
        END DO  
      END DO  
    END DO
    ALLOCATE(A%LUperm(nnzA),LU%LUperm(nnzA))
    A%LUperm  = Permutation
    LU%LUperm = Permutation

  END SUBROUTINE Get_LU_Permutaion
  
  
  
  ! PRINT SPARSE MATRIX (compressed Row format)
  SUBROUTINE PrintSparseMatrix(A)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    !
    INTEGER :: i,j,jj
    !
    !
    WRITE(*,*) 'dimension:     ', A%m,A%n
    WRITE(*,*) 'nonzeros:      ', SIZE(A%ColInd)
    WRITE(*,*)
    ! 
    DO i=1,A%m
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        WRITE(*,'(1X,I5,1X,I5,10X,E23.14)') i,j,A%Val(jj)
      END DO
    END DO
    WRITE(*,*)
  END SUBROUTINE PrintSparseMatrix
  !


  SUBROUTINE WriteSparseMatrix(A,FileName,nr,ns,EXTENDED,CLASSIC,Teq,LinAlg)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    CHARACTER(*), OPTIONAL :: FileName
    INTEGER,      OPTIONAL :: nr,ns
    LOGICAL,      OPTIONAL :: EXTENDED, CLASSIC, Teq
    CHARACTER(2), OPTIONAL :: LinAlg
    !
    INTEGER :: i,j,jj
    !

    ! 
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName))//'_'//TRIM(LinAlg)//'.SparseMat',STATUS='UNKNOWN')

    WRITE(99,*) '###########################################################'
    WRITE(99,*) '##############  Sparse Matrix  Matlab input ###############'
    WRITE(99,*) '###########################################################'
    WRITE(99,*) 
    WRITE(99,*) '###########################################################'
    WRITE(99,*) '#     Name:        ' , ADJUSTL(FileName)
    WRITE(99,*) '#'  
    WRITE(99,*) '#     Dimension:   ' , A%m,' x ',A%n
    WRITE(99,*) '#     Nonzeros:    ' , A%nnz
    WRITE(99,*) '#     nreac, nspc: ' , nr,' , ',ns
    WRITE(99,*) '#     Matrix Form: ' , LinAlg
    WRITE(99,*) '###########################################################'
    WRITE(99,*)

    WRITE(99,*) '###########################################################'
    WRITE(99,*) '# Sparse Matrix in row/column - index format'
    WRITE(99,*) '###########################################################'
    WRITE(99,*) 'MATRIX'
    DO i=1,A%m
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        WRITE(99,'(1X,I12,10X,I12,10X,E23.14)') i,j,A%Val(jj)
      END DO
    END DO
    WRITE(99,*)
    WRITE(99,*)

    IF (ALLOCATED(A%DiagPtr)) THEN
      WRITE(99,*) '###########################################################'
      WRITE(99,*) '# Array for diagonal entries'
      WRITE(99,*) '###########################################################'
      WRITE(99,*) 'DIAG_PTR'
      DO i=1,SIZE(A%DiagPtr)
        WRITE(99,'(1X,I6,10X,I6)') i,A%DiagPtr(i)
      END DO
      WRITE(99,*)
      WRITE(99,*)
    END IF


    IF ( EXTENDED ) THEN
      IF (ALLOCATED(A%DiagPtr_R)) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Array for reaction rate diagonal entries  '
        WRITE(99,*) '#                        (only EXTENDED case) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'DIAG_R_PTR'
        DO i=1,SIZE(A%DiagPtr_R)
          WRITE(99,'(1X,I6,10X,I6)') i,A%DiagPtr_R(i)
        END DO
        WRITE(99,*)
        WRITE(99,*)
      END IF

      IF (ALLOCATED(A%DiagPtr_C)) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Array for species concentration diagonal  '
        WRITE(99,*) '#                        (only EXTENDED case) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'DIAG_C_PTR'
        DO i=1,SIZE(A%DiagPtr_C)
          WRITE(99,'(1X,I6,10X,I6)') i,A%DiagPtr_C(i)
        END DO
        WRITE(99,*)
        WRITE(99,*)
      END IF
    END IF

    IF (ALLOCATED(A%Permu)) THEN
      WRITE(99,*) '##############################################################'
      WRITE(99,*) '# Array for LU permutation and inverse permutation (symmetric)'
      WRITE(99,*) '##############################################################'
      WRITE(99,*) 'DIAG_PERMUTATION'
      DO i=1,SIZE(A%Permu)
        WRITE(99,'(1X,I6,10X,I6,10X,I6)') i,A%Permu(i), A%InvPer(i)
      END DO
      WRITE(99,*)
      WRITE(99,*)
    END IF

    IF (ALLOCATED(A%LUperm)) THEN
      WRITE(99,*) '###########################################################'
      WRITE(99,*) '# Array for LU permutation of all nonzeros '
      WRITE(99,*) '###########################################################'
      WRITE(99,*) 'LU_PERMUTATION'
      DO i=1,SIZE(A%LUperm)
        WRITE(99,'(1X,I6,10X,I6)') i,A%LUperm(i)
      END DO
      WRITE(99,*)
      WRITE(99,*)
    END IF

    IF (Teq) THEN
      IF (ALLOCATED(A%RowVectorPtr)) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Array for indices pointing to the row vector             '
        WRITE(99,*) '#                        (only with Teq) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'ROW_PTR'
        DO i=1,SIZE(A%RowVectorPtr)
          WRITE(99,'(1X,I12,10X,I12)') i,A%RowVectorPtr(i)
        END DO
        WRITE(99,*)
        WRITE(99,*)
      END IF
      IF (ALLOCATED(A%ColVectorPtr)) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Array for indices pointing to the column vector          '
        WRITE(99,*) '#                        (only with Teq) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'COL_PTR'
        DO i=1,SIZE(A%ColVectorPtr)
          WRITE(99,'(1X,I12,10X,I12)') i,A%ColVectorPtr(i)
        END DO
        WRITE(99,*)
        WRITE(99,*)
      END IF
      IF (A%XPtr/=-42) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Index pointing to the X value           '
        WRITE(99,*) '#                        (only with Teq) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'X_PTR'
        WRITE(99,'(1X,I12,10X,I12)') i,A%XPtr
        WRITE(99,*)
        WRITE(99,*)
      ELSE
        WRITE(99,*) 'X_PTR is not definded'
      END IF
    END IF

    CLOSE(99)
    WRITE(*,*) '  Writing matrices to file: ',TRIM(FileName)//'_'//TRIM(LinAlg)//'.SparseMat'
  END SUBROUTINE WriteSparseMatrix


  FUNCTION ReadSparseMatrix(FileName) RESULT(A)
    TYPE(CSR_Matrix_T) :: A
    CHARACTER(*)       :: FileName
    !
    CHARACTER(400) :: io_msg
    INTEGER        :: io_stat
    CHARACTER(400) :: dummy
    INTEGER :: i ,m, n, nnz
    INTEGER,  ALLOCATABLE :: iA(:), jA(:) 
    REAL(dp), ALLOCATABLE :: vA(:)


    ! Synatx:
    !     1. line :: [number of rows] space [number of columns]
    !     2. line :: [number of nonzeros]
    !     3. line :: blank line
    !     4. line - nnz'th line :: [i_row] space [i_column] space [i_val]
    !                             where i=1,...,nnz
    !
    ! 
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName)),STATUS='UNKNOWN',IOSTAT=io_stat,IOMSG=io_msg)

    READ(99,*,IOSTAT=io_stat,IOMSG=io_msg) m, n, dummy
    READ(99,*,IOSTAT=io_stat,IOMSG=io_msg) nnz, dummy
    READ(99,*,IOSTAT=io_stat,IOMSG=io_msg) dummy

    ALLOCATE(iA(nnz),jA(nnz),vA(nnz))
    iA = -99
    jA = -99
    vA = 0.0_dp

    DO i=1,nnz
      READ(99,*,IOSTAT=io_stat,IOMSG=io_msg) iA(i),jA(i),vA(i)
    END DO

    A = New_CSR(m,n,nnz,iA,jA,vA)
    
    CLOSE(99)
    WRITE(*,*)
    WRITE(*,*) '  Read matrix  ::  ',TRIM(FileName), '  done.'
  END FUNCTION ReadSparseMatrix


  FUNCTION InputPivotOrder(FileName,dim) RESULT(PivOrd)
    INTEGER, ALLOCATABLE :: PivOrd(:)
    CHARACTER(*)         :: FileName
    INTEGER              :: dim
    !
    CHARACTER(400) :: io_msg
    INTEGER        :: io_stat, i

    ! Synatx:
    !     1. line - dim'th line :: [i_piv]
    !                             where i=1,...,dim
    !
    ! 
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName)),STATUS='UNKNOWN',IOSTAT=io_stat,IOMSG=io_msg)

    ALLOCATE(PivOrd(dim))
    PivOrd = -99

    DO i=1,dim
      READ(99,*,IOSTAT=io_stat,IOMSG=io_msg) PivOrd(i)
    END DO

    
    CLOSE(99)
    WRITE(*,*)
    WRITE(*,*) '  Read Pivot Order  ::  ',TRIM(FileName), '  done.'
  END FUNCTION InputPivotOrder
  !
  
  ! convert compressed row format to rowindex column index format
  FUNCTION CSR_to_SpRowIndColInd(MatIn) RESULT(MatOut)
    ! only n by n matrices
    TYPE(CSR_Matrix_T) :: MatIn
    TYPE(SpRowIndColInd_T) :: MatOut
    !
    INTEGER :: i,j,jj,n
    MatOut%m=MatIn%m
    MatOut%n=MatIn%n
    !MatOut%n=MatIn%RowPtr(MatIn%m+1)-1

    n = MatIn%RowPtr(MatIn%m+1)-1
    ALLOCATE(MatOut%RowInd(n))
    MatOut%RowInd=0
    ALLOCATE(MatOut%ColInd(n))
    MatOut%ColInd=0
    ALLOCATE(MatOut%Val(n))
    MatOut%Val=0.0_dp
    j=1
    DO i=1,MatIn%m
      DO jj=MatIn%RowPtr(i),MatIn%RowPtr(i+1)-1
        MatOut%RowInd(j)=i
        MatOut%ColInd(j)=MatIn%ColInd(jj)
        MatOut%Val(j)=MatIn%Val(jj)
        j=j+1
      END DO
    END DO

    MatOut%nnz = SIZE(MatOut%RowInd)
  END FUNCTION CSR_to_SpRowIndColInd


  ! sparse matrix * vector + vector
  PURE FUNCTION DAXPY_sparse(A,X,Y) RESULT(Rhs)
    TYPE(CSR_Matrix_T), INTENT(IN)    :: A
    REAL(dp),           INTENT(IN)    :: X(:), Y(:)
    REAL(dp)                          :: Rhs(A%m)
    !TEMP
    INTEGER :: i, rp_i, rp_i1
   
    Rhs = 0.0_dp
    DO i=1,A%m
      rp_i   = A%RowPtr(i)  ;  rp_i1  = A%RowPtr(i+1)-1
      Rhs(i) = SUM(A%Val(rp_i:rp_i1)*X(A%ColInd(rp_i:rp_i1)))+Y(i)
    END DO
  END FUNCTION DAXPY_sparse


  ! sparse matrix * vector
  PURE FUNCTION DAX_sparse(A,X) RESULT(Rhs)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    REAL(dp),           INTENT(IN) :: X(:)
    REAL(dp)                       :: Rhs(A%m)
    !TEMP
    INTEGER :: i, rp_i, rp_i1 
   
    Rhs = 0.0_dp
    DO i=1,A%m
      rp_i   = A%RowPtr(i)
      rp_i1  = A%RowPtr(i+1)-1
      Rhs(i) = SUM(A%Val(rp_i:rp_i1)*X(A%ColInd(rp_i:rp_i1)))
    END DO
  END FUNCTION DAX_sparse
 

  PURE SUBROUTINE SparseLU(A)
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: A
    !
    REAL(dp) :: w(A%n)
    REAL(dp) :: alpha
    INTEGER :: i,j,jj,kk
    !INTEGER :: OPCOUNT

    !OPCOUNT=0
    !
    DO i=1,A%n
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        w(A%ColInd(jj))=A%Val(jj)
      END DO
      DO jj=A%RowPtr(i),A%DiagPtr(i)-1
        j=A%ColInd(jj)
        alpha=w(j)/A%Val(A%DiagPtr(j))
        !OPCOUNT=OPCOUNT+1
        w(j)=alpha
        DO kk=A%DiagPtr(j)+1,A%RowPtr(j+1)-1
          w(A%ColInd(kk))=w(A%ColInd(kk))-alpha*A%Val(kk)
          !OPCOUNT=OPCount+1
        END DO
      END DO
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        A%Val(jj)=w(A%ColInd(jj))
      END DO
    END DO
    !print*, 'opcount=', opcount
    !stop
  END SUBROUTINE SparseLU
  !
  !
  PURE SUBROUTINE SolveSparse(LU,rhs)
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: LU
    REAL(dp),           INTENT(INOUT) :: Rhs(:)
    !
    INTEGER :: i,jj
    REAL(dp) :: b(LU%n)
    
    b( LU%Permu ) = Rhs      
    
    !--  L-solve
    DO i=2,LU%n
      DO jj=LU%RowPtr(i),LU%DiagPtr(i)-1
        b(i)=b(i)-LU%Val(jj)*b(LU%ColInd(jj))
      END DO
    END DO
  
    !--  U-solve
    DO i=LU%n,1,-1
      DO jj=LU%DiagPtr(i)+1,LU%RowPtr(i+1)-1
        b(i)=b(i)-LU%Val(jj)*b(LU%ColInd(jj))
      END DO
      b(i)=b(i)/LU%Val(LU%DiagPtr(i))
    END DO
    
    !--  Back-Permutation of solution
    Rhs( LU%InvPer ) = b
  END SUBROUTINE SolveSparse


END MODULE Sparse_Mod
