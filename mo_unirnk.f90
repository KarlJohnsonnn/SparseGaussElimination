!===============================================================================
! Module: mo_unirnk
! 
! Purpose: Merge-sort ranking of an array with removal of duplicate entries
!
! Description:
!   The routine is similar to pure merge-sort ranking, but on the last pass,
!   it discards indices that correspond to duplicate entries. For performance
!   reasons, the first 2 passes are taken out of the standard loop and use
!   dedicated coding.
!
! Public Interface:
!   unirnk - Generic interface for ranking arrays (double, real, integer)
!===============================================================================
module mo_unirnk

  implicit none
  
  integer, parameter :: kdp = selected_real_kind(15)
  
  public :: unirnk
  private :: kdp
  private :: r_unirnk, i_unirnk, d_unirnk
  private :: r_nearless, i_nearless, d_nearless, nearless
  
  interface unirnk
    module procedure d_unirnk, r_unirnk, i_unirnk
  end interface unirnk
  
  interface nearless
    module procedure d_nearless, r_nearless, i_nearless
  end interface nearless

contains

subroutine d_unirnk(xvalt, irngt, nuni)
  ! Merge-sort ranking for double precision arrays with duplicate removal
  !
  ! Arguments:
  !   xvalt - input array to be ranked
  !   irngt - output permutation vector (rank indices)
  !   nuni  - output number of unique elements
  
  real(kind=8), dimension(:), intent(in) :: xvalt
  integer, dimension(:), intent(out) :: irngt
  integer, intent(out) :: nuni
  
  integer, dimension(size(irngt)) :: jwrkt
  integer :: lmtna, lmtnc, irng, irng1, irng2
  integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
  real(kind=8) :: xtst, xvala, xvalb

  nval = min(size(xvalt), size(irngt))
  nuni = nval

  select case (nval)
  case (:0)
    return
  case (1)
    irngt(1) = 1
    return
  end select

  ! Fill-in the index array, creating ordered couples
  do iind = 2, nval, 2
    if (xvalt(iind-1) < xvalt(iind)) then
      irngt(iind-1) = iind - 1
      irngt(iind) = iind
    else
      irngt(iind-1) = iind
      irngt(iind) = iind - 1
    end if
  end do
  
  if (modulo(nval, 2) /= 0) then
    irngt(nval) = nval
  end if

  lmtna = 2
  lmtnc = 4

  ! First iteration: ordered subset length goes from 2 to 4
  do
    if (nval <= 4) exit

    do iwrkd = 0, nval - 1, 4
      if ((iwrkd+4) > nval) then
        if ((iwrkd+2) >= nval) exit
        if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit

        if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
          irng2 = irngt(iwrkd+2)
          irngt(iwrkd+2) = irngt(iwrkd+3)
          irngt(iwrkd+3) = irng2
        else
          irng1 = irngt(iwrkd+1)
          irngt(iwrkd+1) = irngt(iwrkd+3)
          irngt(iwrkd+3) = irngt(iwrkd+2)
          irngt(iwrkd+2) = irng1
        end if
        exit
      end if

      if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle

      if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
        irng2 = irngt(iwrkd+2)
        irngt(iwrkd+2) = irngt(iwrkd+3)
        if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
          irngt(iwrkd+3) = irng2
        else
          irngt(iwrkd+3) = irngt(iwrkd+4)
          irngt(iwrkd+4) = irng2
        end if
      else
        irng1 = irngt(iwrkd+1)
        irng2 = irngt(iwrkd+2)
        irngt(iwrkd+1) = irngt(iwrkd+3)
        if (xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
          irngt(iwrkd+2) = irng1
          if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
            irngt(iwrkd+3) = irng2
          else
            irngt(iwrkd+3) = irngt(iwrkd+4)
            irngt(iwrkd+4) = irng2
          end if
        else
          irngt(iwrkd+2) = irngt(iwrkd+4)
          irngt(iwrkd+3) = irng1
          irngt(iwrkd+4) = irng2
        end if
      end if
    end do

    lmtna = 4
    exit
  end do

  ! Main iteration loop: double the length of ordered subsets each time
  do
    if (2*lmtna >= nval) exit
    iwrkf = 0
    lmtnc = 2 * lmtnc

    do
      iwrk = iwrkf
      iwrkd = iwrkf + 1
      jinda = iwrkf + lmtna
      iwrkf = iwrkf + lmtnc
      if (iwrkf >= nval) then
        if (jinda >= nval) exit
        iwrkf = nval
      end if
      iinda = 1
      iindb = jinda + 1

      jwrkt(1:lmtna) = irngt(iwrkd:jinda)
      xvala = xvalt(jwrkt(iinda))
      xvalb = xvalt(irngt(iindb))

      do
        iwrk = iwrk + 1

        if (xvala > xvalb) then
          irngt(iwrk) = irngt(iindb)
          iindb = iindb + 1
          if (iindb > iwrkf) then
            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
            exit
          end if
          xvalb = xvalt(irngt(iindb))
        else
          irngt(iwrk) = jwrkt(iinda)
          iinda = iinda + 1
          if (iinda > lmtna) exit
          xvala = xvalt(jwrkt(iinda))
        end if
      end do
    end do

    lmtna = 2 * lmtna
  end do

  ! Last merge with removal of duplicates
  iinda = 1
  iindb = lmtna + 1
  nuni = 0

  jwrkt(1:lmtna) = irngt(1:lmtna)
  if (iindb <= nval) then
    xtst = nearless(min(xvalt(jwrkt(1)), xvalt(irngt(iindb))))
  else
    xtst = nearless(xvalt(jwrkt(1)))
  end if
  
  do iwrk = 1, nval
    if (iinda <= lmtna) then
      if (iindb <= nval) then
        if (xvalt(jwrkt(iinda)) > xvalt(irngt(iindb))) then
          irng = irngt(iindb)
          iindb = iindb + 1
        else
          irng = jwrkt(iinda)
          iinda = iinda + 1
        end if
      else
        irng = jwrkt(iinda)
        iinda = iinda + 1
      end if
    else
      irng = irngt(iwrk)
    end if
    
    if (xvalt(irng) > xtst) then
      xtst = xvalt(irng)
      nuni = nuni + 1
      irngt(nuni) = irng
    end if
  end do

end subroutine d_unirnk

subroutine r_unirnk(xvalt, irngt, nuni)
  ! Merge-sort ranking for single precision arrays with duplicate removal
  
  real(kind=4), dimension(:), intent(in) :: xvalt
  integer, dimension(:), intent(out) :: irngt
  integer, intent(out) :: nuni
  
  integer, dimension(size(irngt)) :: jwrkt
  integer :: lmtna, lmtnc, irng, irng1, irng2
  integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
  real(kind=4) :: xtst, xvala, xvalb

  nval = min(size(xvalt), size(irngt))
  nuni = nval

  select case (nval)
  case (:0)
    return
  case (1)
    irngt(1) = 1
    return
  end select

  do iind = 2, nval, 2
    if (xvalt(iind-1) < xvalt(iind)) then
      irngt(iind-1) = iind - 1
      irngt(iind) = iind
    else
      irngt(iind-1) = iind
      irngt(iind) = iind - 1
    end if
  end do
  
  if (modulo(nval, 2) /= 0) then
    irngt(nval) = nval
  end if

  lmtna = 2
  lmtnc = 4

  do
    if (nval <= 4) exit

    do iwrkd = 0, nval - 1, 4
      if ((iwrkd+4) > nval) then
        if ((iwrkd+2) >= nval) exit
        if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit

        if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
          irng2 = irngt(iwrkd+2)
          irngt(iwrkd+2) = irngt(iwrkd+3)
          irngt(iwrkd+3) = irng2
        else
          irng1 = irngt(iwrkd+1)
          irngt(iwrkd+1) = irngt(iwrkd+3)
          irngt(iwrkd+3) = irngt(iwrkd+2)
          irngt(iwrkd+2) = irng1
        end if
        exit
      end if

      if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle

      if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
        irng2 = irngt(iwrkd+2)
        irngt(iwrkd+2) = irngt(iwrkd+3)
        if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
          irngt(iwrkd+3) = irng2
        else
          irngt(iwrkd+3) = irngt(iwrkd+4)
          irngt(iwrkd+4) = irng2
        end if
      else
        irng1 = irngt(iwrkd+1)
        irng2 = irngt(iwrkd+2)
        irngt(iwrkd+1) = irngt(iwrkd+3)
        if (xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
          irngt(iwrkd+2) = irng1
          if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
            irngt(iwrkd+3) = irng2
          else
            irngt(iwrkd+3) = irngt(iwrkd+4)
            irngt(iwrkd+4) = irng2
          end if
        else
          irngt(iwrkd+2) = irngt(iwrkd+4)
          irngt(iwrkd+3) = irng1
          irngt(iwrkd+4) = irng2
        end if
      end if
    end do

    lmtna = 4
    exit
  end do

  do
    if (2*lmtna >= nval) exit
    iwrkf = 0
    lmtnc = 2 * lmtnc

    do
      iwrk = iwrkf
      iwrkd = iwrkf + 1
      jinda = iwrkf + lmtna
      iwrkf = iwrkf + lmtnc
      if (iwrkf >= nval) then
        if (jinda >= nval) exit
        iwrkf = nval
      end if
      iinda = 1
      iindb = jinda + 1

      jwrkt(1:lmtna) = irngt(iwrkd:jinda)
      xvala = xvalt(jwrkt(iinda))
      xvalb = xvalt(irngt(iindb))

      do
        iwrk = iwrk + 1

        if (xvala > xvalb) then
          irngt(iwrk) = irngt(iindb)
          iindb = iindb + 1
          if (iindb > iwrkf) then
            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
            exit
          end if
          xvalb = xvalt(irngt(iindb))
        else
          irngt(iwrk) = jwrkt(iinda)
          iinda = iinda + 1
          if (iinda > lmtna) exit
          xvala = xvalt(jwrkt(iinda))
        end if
      end do
    end do

    lmtna = 2 * lmtna
  end do

  iinda = 1
  iindb = lmtna + 1
  nuni = 0

  jwrkt(1:lmtna) = irngt(1:lmtna)
  if (iindb <= nval) then
    xtst = nearless(min(xvalt(jwrkt(1)), xvalt(irngt(iindb))))
  else
    xtst = nearless(xvalt(jwrkt(1)))
  end if
  
  do iwrk = 1, nval
    if (iinda <= lmtna) then
      if (iindb <= nval) then
        if (xvalt(jwrkt(iinda)) > xvalt(irngt(iindb))) then
          irng = irngt(iindb)
          iindb = iindb + 1
        else
          irng = jwrkt(iinda)
          iinda = iinda + 1
        end if
      else
        irng = jwrkt(iinda)
        iinda = iinda + 1
      end if
    else
      irng = irngt(iwrk)
    end if
    
    if (xvalt(irng) > xtst) then
      xtst = xvalt(irng)
      nuni = nuni + 1
      irngt(nuni) = irng
    end if
  end do

end subroutine r_unirnk

subroutine i_unirnk(xvalt, irngt, nuni)
  ! Merge-sort ranking for integer arrays with duplicate removal
  
  integer, dimension(:), intent(in) :: xvalt
  integer, dimension(:), intent(out) :: irngt
  integer, intent(out) :: nuni
  
  integer, dimension(size(irngt)) :: jwrkt
  integer :: lmtna, lmtnc, irng, irng1, irng2
  integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
  integer :: xtst, xvala, xvalb

  nval = min(size(xvalt), size(irngt))
  nuni = nval

  select case (nval)
  case (:0)
    return
  case (1)
    irngt(1) = 1
    return
  end select

  do iind = 2, nval, 2
    if (xvalt(iind-1) < xvalt(iind)) then
      irngt(iind-1) = iind - 1
      irngt(iind) = iind
    else
      irngt(iind-1) = iind
      irngt(iind) = iind - 1
    end if
  end do
  
  if (modulo(nval, 2) /= 0) then
    irngt(nval) = nval
  end if

  lmtna = 2
  lmtnc = 4

  do
    if (nval <= 4) exit

    do iwrkd = 0, nval - 1, 4
      if ((iwrkd+4) > nval) then
        if ((iwrkd+2) >= nval) exit
        if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) exit

        if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
          irng2 = irngt(iwrkd+2)
          irngt(iwrkd+2) = irngt(iwrkd+3)
          irngt(iwrkd+3) = irng2
        else
          irng1 = irngt(iwrkd+1)
          irngt(iwrkd+1) = irngt(iwrkd+3)
          irngt(iwrkd+3) = irngt(iwrkd+2)
          irngt(iwrkd+2) = irng1
        end if
        exit
      end if

      if (xvalt(irngt(iwrkd+2)) <= xvalt(irngt(iwrkd+3))) cycle

      if (xvalt(irngt(iwrkd+1)) <= xvalt(irngt(iwrkd+3))) then
        irng2 = irngt(iwrkd+2)
        irngt(iwrkd+2) = irngt(iwrkd+3)
        if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
          irngt(iwrkd+3) = irng2
        else
          irngt(iwrkd+3) = irngt(iwrkd+4)
          irngt(iwrkd+4) = irng2
        end if
      else
        irng1 = irngt(iwrkd+1)
        irng2 = irngt(iwrkd+2)
        irngt(iwrkd+1) = irngt(iwrkd+3)
        if (xvalt(irng1) <= xvalt(irngt(iwrkd+4))) then
          irngt(iwrkd+2) = irng1
          if (xvalt(irng2) <= xvalt(irngt(iwrkd+4))) then
            irngt(iwrkd+3) = irng2
          else
            irngt(iwrkd+3) = irngt(iwrkd+4)
            irngt(iwrkd+4) = irng2
          end if
        else
          irngt(iwrkd+2) = irngt(iwrkd+4)
          irngt(iwrkd+3) = irng1
          irngt(iwrkd+4) = irng2
        end if
      end if
    end do

    lmtna = 4
    exit
  end do

  do
    if (2*lmtna >= nval) exit
    iwrkf = 0
    lmtnc = 2 * lmtnc

    do
      iwrk = iwrkf
      iwrkd = iwrkf + 1
      jinda = iwrkf + lmtna
      iwrkf = iwrkf + lmtnc
      if (iwrkf >= nval) then
        if (jinda >= nval) exit
        iwrkf = nval
      end if
      iinda = 1
      iindb = jinda + 1

      jwrkt(1:lmtna) = irngt(iwrkd:jinda)
      xvala = xvalt(jwrkt(iinda))
      xvalb = xvalt(irngt(iindb))

      do
        iwrk = iwrk + 1

        if (xvala > xvalb) then
          irngt(iwrk) = irngt(iindb)
          iindb = iindb + 1
          if (iindb > iwrkf) then
            irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
            exit
          end if
          xvalb = xvalt(irngt(iindb))
        else
          irngt(iwrk) = jwrkt(iinda)
          iinda = iinda + 1
          if (iinda > lmtna) exit
          xvala = xvalt(jwrkt(iinda))
        end if
      end do
    end do

    lmtna = 2 * lmtna
  end do

  iinda = 1
  iindb = lmtna + 1
  nuni = 0

  jwrkt(1:lmtna) = irngt(1:lmtna)
  if (iindb <= nval) then
    xtst = nearless(min(xvalt(jwrkt(1)), xvalt(irngt(iindb))))
  else
    xtst = nearless(xvalt(jwrkt(1)))
  end if
  
  do iwrk = 1, nval
    if (iinda <= lmtna) then
      if (iindb <= nval) then
        if (xvalt(jwrkt(iinda)) > xvalt(irngt(iindb))) then
          irng = irngt(iindb)
          iindb = iindb + 1
        else
          irng = jwrkt(iinda)
          iinda = iinda + 1
        end if
      else
        irng = jwrkt(iinda)
        iinda = iinda + 1
      end if
    else
      irng = irngt(iwrk)
    end if
    
    if (xvalt(irng) > xtst) then
      xtst = xvalt(irng)
      nuni = nuni + 1
      irngt(nuni) = irng
    end if
  end do

end subroutine i_unirnk

function d_nearless(xval) result(d_nl)
  ! Return nearest value less than given double precision value
  real(kind=8), intent(in) :: xval
  real(kind=8) :: d_nl
  
  d_nl = nearest(xval, -1.0_kdp)
  
end function d_nearless

function r_nearless(xval) result(r_nl)
  ! Return nearest value less than given single precision value
  real(kind=4), intent(in) :: xval
  real(kind=4) :: r_nl
  
  r_nl = nearest(xval, -1.0)
  
end function r_nearless

function i_nearless(xval) result(i_nl)
  ! Return nearest value less than given integer value
  integer, intent(in) :: xval
  integer :: i_nl
  
  i_nl = xval - 1
  
end function i_nearless

end module mo_unirnk
