module jack

  use types, only: WP

  implicit none

  public :: Complement

  interface Complement
     module procedure Complement0
     module procedure Complement1
     module procedure Complement2
     module procedure Complement3
  end interface Complement

contains

  pure subroutine Complement0 ( cset, data )
    real(kind=WP),               intent(out) :: cset
    real(kind=WP), dimension(:), intent(in)  :: data
    cset = sum(data) / real(size(data), kind=WP)
  end subroutine Complement0

  pure subroutine Complement1 (ncon, cset, data )
    integer, intent(in) :: ncon
    real(kind=WP), dimension(ncon), intent(out) :: cset
    real(kind=WP), dimension(ncon), intent(in)  :: data
    cset = (sum(data) - data) / real(size(data)-1, kind=WP)
  end subroutine Complement1

  pure subroutine Complement2 ( ncon, cset, data, icon )
    integer, intent(in) :: ncon
    real(kind=WP), dimension(ncon), intent(out) :: cset
    real(kind=WP), dimension(ncon), intent(in)  :: data
    integer,                     intent(in)  :: icon
    cset = (sum(data) - data(icon) - data) / real(size(data)-2, kind=WP)
  end subroutine Complement2

  pure subroutine Complement3 ( ncon, cset, data, icon, jcon )
    integer, intent(in) :: ncon
    real(kind=WP), dimension(ncon), intent(out) :: cset
    real(kind=WP), dimension(ncon), intent(in)  :: data
    integer,                     intent(in)  :: icon
    integer,                     intent(in)  :: jcon
    cset = (sum(data) - data(icon) - data(jcon) - data) / real(size(data)-3, kind=WP)
  end subroutine Complement3

  subroutine JackKnife_wp ( ncon, c, err, bias )
    ! This subroutine works
    ! but just rescale np.cov by (ncon - 1)
    integer, intent(in) :: ncon
    real(kind=WP), dimension(0:ncon), intent(in)  :: c
    real(kind=WP),                intent(out) :: err
    real(kind=WP), optional,      intent(out) :: bias
    real(kind=WP)                             :: avg
    avg = sum(c(1:))/real(size(c)-1, kind=WP)
    err = sqrt(sum((c(1:) - avg)**2)*real(size(c)-2, kind=WP)/real(size(c)-1, kind=WP))
    if (present(bias)) bias = c(0) - avg
  end subroutine JackKnife_wp


end module jack
