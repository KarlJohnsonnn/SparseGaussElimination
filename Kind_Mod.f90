!===============================================================================
! Module: kind_mod
! 
! Purpose: Define precision parameters for floating-point arithmetic
!
! Author: Willi Schimmel
! Institute: Leibniz Institute for Tropospheric Research (TROPOS)
!===============================================================================
module kind_mod

  implicit none
  
  ! Double precision kind parameter
  integer, parameter :: dp = selected_real_kind(15, 307)

end module kind_mod
