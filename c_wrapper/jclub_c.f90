module jclub_c
  use jclub, only: c_double_complex, c_int, c_double, c_char, plaquette, polyakov, &
  readgaugefieldunopenqcd, poly_kolaorder, plaquette_kolaorder, &
  readgaugefieldunopenqcdunkolaorder, complement, spacingCalcs, w0Calcs, xigCalcs
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


  ! Jack
  pure subroutine complement0_c(cset, data)
    real(kind=c_double), intent(out) :: cset
    real(kind=c_double), dimension(:), intent(in) :: data
    call complement(cset, data)
  end subroutine complement0_c

  pure subroutine complement1_c(ncon, cset, data )
    integer(kind=c_int), intent(in) :: ncon
    real(kind=c_double), dimension(ncon), intent(out) :: cset
    real(kind=c_double), dimension(ncon), intent(in)  :: data
    call complement(ncon, cset, data)
  end subroutine complement1_c

  ! flowAna

  subroutine spacingCalcs_c(ncon, xiNum, flowTimeForW0, xiNumList, xig, nEval, &
       w0PhysMean, w0PhysErr, a_sSplineEval, xEval, a_s, a_sSys)
    integer(kind=c_int), intent(in) :: ncon, xiNum
    real(kind=c_double), dimension(xiNum, 0:ncon), intent(in) :: flowTimeForW0 ! xi, 0:ncon
    real(kind=c_double), dimension(xiNum), intent(in) :: xiNumList ! xi
    real(kind=c_double), dimension(0:ncon), intent(in) :: xig ! 0:ncon
    integer(kind=c_int), intent(in) :: nEval
    real(kind=c_double), intent(in) :: w0PhysMean, w0PhysErr  ! in fm
    real(kind=c_double), dimension(nEval, 0:ncon), intent(out) :: a_sSplineEval  ! nEval, 0:ncon
    real(kind=c_double), dimension(nEval), intent(out) :: xEval ! nEval
    real(kind=c_double), dimension(0:ncon), intent(out) :: a_s
    real(kind=c_double), intent(out) :: a_sSys
    call spacingCalcs(flowTimeForW0, xiNumList, xig, nEval, w0PhysMean, w0PhysErr, a_sSplineEval, xEval, a_s, a_sSys)
  end subroutine spacingCalcs_c

  subroutine w0Calcs_c(ncon, xiNum, flowTime, JE4i, JEij, xiNumList, &
       nEval, targws, xEval, wijSplineEval, w4iSplineEval, flowTimeForW0, w0ij, w04i)
    integer(kind=c_int), intent(in) :: ncon, xiNum
    real(kind=c_double), dimension(nEval), intent(in) :: flowTime
    real(kind=c_double), dimension(nEval,xiNum,0:ncon), intent(in) :: JE4i, JEij  ! t, xi, 0:ncon
    real(kind=c_double), dimension(xiNum), intent(in) :: xiNumList
    integer(kind=c_int), intent(in) :: nEval
    real(kind=c_double), intent(in) :: targws
    real(kind=c_double), dimension(nEval,xiNum,0:ncon),  intent(out) :: wijSplineEval, w4iSplineEval ! t, xi, 0:ncon
    real(kind=c_double), dimension(xiNum, 0:ncon), intent(out) :: flowTimeForW0 ! xi, 0:ncon
    real(kind=c_double), dimension(xiNum,0:ncon), intent(out) :: w0ij, w04i ! xi, 0:ncon
    real(kind=c_double), dimension(nEval), intent(out) :: xEval  ! t

    call w0Calcs(flowTime, JE4i, JEij, xiNumList, nEval, targws, xEval, wijSplineEval, w4iSplineEval, flowTimeForW0, w0ij, w04i)
  end subroutine w0Calcs_c

  subroutine xigCalcs_c(xiNum, ncon, RE, xiNumList, nEval, xEval, RESplineEval, xig)
    integer(kind=c_int), intent(in) :: xiNum, ncon
    real(kind=c_double), dimension(xiNum, 0:ncon), intent(in) :: RE  ! xi, 0:ncon
    real(kind=c_double), dimension(xiNum), intent(in) :: xiNumList
    integer(kind=c_int), intent(in) :: neVal
    real(kind=c_double), dimension(nEval), intent(out) :: xEval
    real(kind=c_double), dimension(nEval, 0:ncon), intent(out) :: ReSplineEval ! nEval, 0:ncon
    real(kind=c_double), dimension(0:ncon), intent(out) :: xig ! 0:ncon
    call xigCalcs(RE, xiNumList, nEval, xEval, RESplineEval, xig)
  end subroutine xigCalcs_c

end module jclub_c
