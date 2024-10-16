module jclub_c

  use jclub, only: c_double_complex, c_int, c_double, c_char, plaquette, polyakov, &
  readgaugefieldunopenqcd, poly_kolaorder, plaquette_kolaorder, &
  readgaugefieldunopenqcdunkolaorder
  implicit none

contains

  ! wloops
  subroutine plaquette_c(data, mustart, muend, nuend, sumtrp, np, time)
    complex(kind=c_double_complex), dimension(:, :, :, :, :, :, :), intent(in)  :: data
    integer(kind=c_int),                                          intent(in)  :: mustart, muend, nuend
    real(kind=c_double),                                    intent(out) :: sumtrp, time
    integer(kind=c_int),                                          intent(out) :: np
    call plaquette(data, mustart, muend, nuend, sumtrp, np, time)
  end subroutine plaquette_c

  subroutine polyakov_c(data, sumtrp, np, time)
    complex(kind=c_double_complex), dimension(:, :, :, :, :, :, :), intent(in)  :: data
    real(kind=c_double),                                    intent(out) :: sumtrp, time
    integer(kind=c_int),                                          intent(out) :: np
    call polyakov(data, sumtrp, np, time)
  end subroutine polyakov_c

  subroutine poly_kolaorder_c(data, sumtrp, np, time)
    complex(kind=c_double_complex), dimension(:, :, :, :, :, :, :), intent(in)  :: data
    real(kind=c_double),                                    intent(out) :: sumtrp, time
    integer(kind=c_int),                                          intent(out) :: np
    integer(kind=c_int) :: ns, nt
    integer(kind=c_int), dimension(7) :: datashape
    datashape = shape(data)
    ns = datashape(5)
    nt = datashape(7)
    call poly_kolaorder(data, ns, nt, sumtrp, np, time)
  end subroutine poly_kolaorder_c

  subroutine plaquette_kolaorder_c(data, mustart, muend, nuend, sumtrp, np, time)
    complex(kind=c_double_complex), dimension(:, :, :, :, :, :, :), intent(in)  :: data
    integer(kind=c_int),                                          intent(in)  :: mustart, muend, nuend
    real(kind=c_double),                                    intent(out) :: sumtrp, time
    integer(kind=c_int),                                          intent(out) :: np
    call plaquette_kolaorder(data, mustart, muend, nuend, sumtrp, np, time)
  end subroutine plaquette_kolaorder_c

  ! openqcdio
  function readgaugefieldUNopenqcd_c(filename, nx, ny, nz, nt) result(u_xd)
    character(len=*,kind=c_char), intent(in) :: filename
    integer(kind=c_int),          intent(in) :: nx, ny, nz, nt
    complex(kind=c_double_complex), dimension(nt,nx,ny,nz,4,3,3) :: u_xd
    u_xd = readgaugefieldUNopenqcd(filename, nx, ny, nz, nt)
  end function readgaugefieldUNopenqcd_c

  function readgaugefieldUNopenqcdUNkolaorder_c(filename, nx, ny, nz, nt) result(u_xd)
    character(len=*,kind=c_char), intent(in) :: filename
    integer(kind=c_int),          intent(in) :: nx, ny, nz, nt
    complex(kind=c_double_complex), dimension(3,3,4,nx,nx,ny,nt) :: u_xd
    u_xd = readgaugefieldUNopenqcdUNkolaorder(filename, nx, ny, nz, nt)
  end function readgaugefieldUNopenqcdUNkolaorder_c


end module jclub_c
