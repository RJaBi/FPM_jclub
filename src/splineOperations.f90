module splineOperations
  use types, only: WP
  use bspline_module
  !use root_module
  implicit none


  private

  public :: setupSplineParams, finaliseSplineParams
  public :: ff, minFunc, splineXDerivMinFunc, setSplineTarget
  public :: splineMinFunc

  public :: nx, kx, nxv, bcoef, tx, fval, w1_1d, extrap, xval

  ! Spline parameters...
  integer :: nx, kx, nxv  ! kx = 4 is cubic spline
  real(kind=WP), dimension(:), allocatable :: x, xval, fval, fdval, tx, bcoef, w1_1d
  logical, parameter :: extrap=.true.
  real(kind=WP) :: splineTarg  ! splineTarg is subtracted off

  abstract interface
     function fit_function(x)
       import :: WP
       implicit none
       real(kind=WP), intent(in)    :: x
       real(kind=WP) :: fit_function
     end function fit_function
  end interface

  procedure(fit_function), pointer :: ff

contains

  function minFunc(x) result(f)
    real(kind=WP), intent(in) :: x
    real(kind=WP) :: f
    f = ff(x)
  end function minFunc


  function splineMinFunc(x) result(f)
    real(kind=WP), intent(in) :: x
    real(kind=WP) :: f
    integer :: iflag, inbvx
    call db1val(x, 0, tx, nx, kx, bcoef, f, iflag, inbvx, w1_1d, extrap=extrap)
    f = f - splineTarg
  end function splineMinFunc


  function splineXDerivMinFunc(x) result(f)
    real(kind=WP), intent(in) :: x
    real(kind=WP) :: f
    integer :: iflag, inbvx
    call db1val(x, 1, tx, nx, kx, bcoef, f, iflag, inbvx, w1_1d, extrap=extrap)
    f = x * f - splineTarg
    !f = db1val(x, 1, tx, nx, kx, bcoef, val, iflag, inbvx, w1_1d, extrap=extrap)
  end function splineXDerivMinFunc

  subroutine setSplineTarget(targ)
    real(kind=WP), intent(in) :: targ
    splineTarg = targ
  end subroutine setSplineTarget

  subroutine setupSplineParams(nix, kix, nixv)
    integer, intent(in) :: nix, kix, nixv
    nx = nix
    kx = kix
    nxv = nixv
    allocate(x(nx))
    allocate(xval(nxv))
    allocate(fval(nxv))
    allocate(tx(nx+kx))
    allocate(bcoef(nx))
    allocate(w1_1d(3*kx))
  end subroutine setupSplineParams

  subroutine finaliseSplineParams()
    deallocate(x, xval, fval, tx, bcoef, w1_1d)
  end subroutine finaliseSplineParams


end module splineOperations
