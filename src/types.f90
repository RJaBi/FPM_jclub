module types
  use ISO_C_BINDING, only: C_DOUBLE, C_INT, C_DOUBLE_COMPLEX, C_BOOL, C_CHAR
  implicit none
  ! private
  !! integer, parameter :: dc = kind((1.0D0,1.0D0)) !! Double precision complex scalars
  !integer, parameter :: dp = kind(1.0D0) !! Double precision real
  integer, parameter :: dp = C_DOUBLE
  integer, parameter :: dc = C_DOUBLE_COMPLEX


  integer, parameter, public :: WP = DP
  integer, parameter, public :: WC = DC

  ! real(dc) :: r_dc = (1.0_dp, 1.0_dp)
  ! real(dp) :: r_dp = 1.0_dp

  type :: raggedIntArr
     integer, dimension(:), allocatable :: rag
  end type raggedIntArr

end module types
