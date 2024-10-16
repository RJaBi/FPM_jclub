!! Calculate, spatial, temporal plaquettes
!! Compiles for use with python using
!! 'f2py -c fortPlaq.f90 -m fortPlaq'
!!
!! Which can then be imported into python as
!! import fortPlaq
!! fPlaq = fortPlaq.plaq.plaquette
!! which would give you the plaquette function



module wloops
  use types, only: WP, WC , c_bool
  !use openQCDIO
  ! use omp_lib
  implicit none

  !private
  public

  !public :: plaquette, polyakov

  !public :: plaquette_kolaorder, poly_kolaorder

contains

  pure subroutine MultiplyMatMat(MM, left, right)
    complex(kind=WC), dimension(3, 3), intent(in) :: left, right
    complex(kind=WC), dimension(3, 3), intent(out) :: MM
    !"""
    !Multiple left by right. Assumes 3x3 (colour) (complex) matrices
    !"""
    !# do the maths for the colour matrices
    MM(1, 1) = left(1, 1) * right(1, 1) + left(1, 2) * right(2, 1) + left(1, 3) * right(3, 1)
    MM(2, 1) = left(2, 1) * right(1, 1) + left(2, 2) * right(2, 1) + left(2, 3) * right(3, 1)
    MM(3, 1) = left(3, 1) * right(1, 1) + left(3, 2) * right(2, 1) + left(3, 3) * right(3, 1)
    !# second index
    MM(1, 2) = left(1, 1) * right(1, 2) + left(1, 2) * right(2, 2) + left(1, 3) * right(3, 2)
    MM(2, 2) = left(2, 1) * right(1, 2) + left(2, 2) * right(2, 2) + left(2, 3) * right(3, 2)
    MM(3, 2) = left(3, 1) * right(1, 2) + left(3, 2) * right(2, 2) + left(3, 3) * right(3, 2)
    !# third index
    MM(1, 3) = left(1, 1) * right(1, 3) + left(1, 2) * right(2, 3) + left(1, 3) * right(3, 3)
    MM(2, 3) = left(2, 1) * right(1, 3) + left(2, 2) * right(2, 3) + left(2, 3) * right(3, 3)
    MM(3, 3) = left(3, 1) * right(1, 3) + left(3, 2) * right(2, 3) + left(3, 3) * right(3, 3)
  end subroutine MultiplyMatMat


  pure subroutine MultiplyMatdagMatdag(MM, left, right)
    complex(kind=WC), dimension(3, 3), intent(in) :: left, right
    complex(kind=WC), dimension(3, 3), intent(out) :: MM
    !"""
    !#Multiplies two (3,3) complex matrices together. Takes conjugate
    !Does (left*right)^dagger
    !"""
    !# take transpose manually
    MM(1, 1) = conjg(left(1, 1) * right(1, 1) + left(2, 1) * right(1, 2) + left(3, 1) * right(1, 3))
    MM(2, 1) = conjg(left(1, 2) * right(1, 1) + left(2, 2) * right(1, 2) + left(3, 2) * right(1, 3))
    MM(3, 1) = conjg(left(1, 3) * right(1, 1) + left(2, 3) * right(1, 2) + left(3, 3) * right(1, 3))
    !# but take conjugate using np
    MM(1, 2) = conjg(left(1, 1) * right(2, 1) + left(2, 1) * right(2, 2) + left(3, 1) * right(2, 3))
    MM(2, 2) = conjg(left(1, 2) * right(2, 1) + left(2, 2) * right(2, 2) + left(3, 2) * right(2, 3))
    MM(3, 2) = conjg(left(1, 3) * right(2, 1) + left(2, 3) * right(2, 2) + left(3, 3) * right(2, 3))
    !# last index
    MM(1, 3) = conjg(left(1, 1) * right(3, 1) + left(2, 1) * right(3, 2) + left(3, 1) * right(3, 3))
    MM(2, 3) = conjg(left(1, 2) * right(3, 1) + left(2, 2) * right(3, 2) + left(3, 2) * right(3, 3))
    MM(3, 3) = conjg(left(1, 3) * right(3, 1) + left(2, 3) * right(3, 2) + left(3, 3) * right(3, 3))
  end subroutine MultiplyMatdagMatdag


  pure subroutine RealTraceMultMatMat(TrMM, left, right)
    complex(kind=WC), dimension(3, 3), intent(in) :: left, right
    real(kind=WC), intent(out) :: TrMM
    !"""
    !# !Takes the real trace of (3,3) complex numbers left, right multiplied together
    !Tr(left*right)
    !"""
    TrMM = real(left(1, 1) * right(1, 1) + left(1, 2) * right(2, 1) + left(1, 3) * right(3, 1) + &
         left(2, 1) * right(1, 2) + left(2, 2) * right(2, 2) + left(2, 3) * right(3, 2) + &
         left(3, 1) * right(1, 3) + left(3, 2) * right(2, 3) + left(3, 3) * right(3, 3), kind=WP)
  end subroutine RealTraceMultMatMat

  pure function set3x3Ident()
    complex (kind=WC), dimension(3,3) ::set3x3Ident
    ! A short function to constrcuct the compelx 3x3 identity matrix
    set3x3Ident = (0.0, 0.0)
    set3x3Ident(1, 1) = (1.0, 0.0)
    set3x3Ident(2, 2) = (1.0, 0.0)
    set3x3Ident(3, 3) = (1.0, 0.0)
  end function set3x3Ident

  subroutine plaquette(data, muStart, muEnd, nuEnd, sumTrP, nP, time)
    complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in)  :: data
    integer,                                          intent(in)  :: muStart, muEnd, nuEnd
    real(kind=WP),                                    intent(out) :: sumTrP, time
    integer,                                          intent(out) :: nP
    !"""
    !Calculates the plaquette over muStart to muEnd
    !data is [nt, nx, ny, nz, mu, colour, colour] complex
    !the plaquette over all lattice is muStart=1, muEnd=4, nuEnd=4
    !the spatial plaquette is muStart=2, muEnd=4, nuEnd=4
    !the temporal plaquette is muStart=1, muEnd=1, nuEnd=4
    !returns the sum of plaquettes, number of plaquettes measured,
    !the average plaquette and the time taken to calculate it
    !"""
    integer, dimension(7) :: dataShape
    integer, dimension(4) :: muCoord, nuCoord, coordBase, coord
    integer :: mu, nu, nx, ny, nz, nt, cc  ! Counters
    ! For intermediate calculating plaquette
    complex(kind=WC), dimension(3, 3) :: Umu_x, Unu_xmu, UmuUnu
    complex(kind=WC), dimension(3, 3) :: Umu_xnu, Unu_x, UmudagUnudag
    real(kind=WP) :: P
    ! Timers
    real :: start, end
    call cpu_time(start)
    dataShape = shape(data)
    !# hold the sum
    sumTrP = 0.0
    !# hold the number measured
    nP = 0
    do mu = muStart, muEnd
       muCoord(:) = 0
       !# This is the shift in mu
       muCoord(mu) = 1
       do nu = mu + 1, nuEnd
          nuCoord(:) = 0
          !# This is the shift in nu
          nuCoord(nu) = 1
          !# loop over all sites
          do nx = 1, dataShape(2)
             do ny = 1, dataShape(3)
                do nz = 1, dataShape(4)
                   do nt = 1, dataShape(1)
                      !# U_mu(x)
                      coordBase = (/ nt, nx, ny, nz /)
                      coord = coordBase
                      Umu_x = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                      !# U_nu(x + amu)
                      coord = coordBase + muCoord
                      !# respect periodic boundary conditions
                      do cc = 1, size(coord)
                         if (coord(cc) > dataShape(cc)) then
                            coord(cc) = 1
                         end if
                      end do
                      Unu_xmu = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                      !# U_mu(x + anu)
                      coord = coordBase + nuCoord
                      do cc = 1, size(coord)
                         if (coord(cc) > dataShape(cc)) then
                            coord(cc) = 1
                         end if
                      end do
                      Umu_xnu = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                      !# U_nu(x)
                      coord = coordBase
                      Unu_x = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                      !# Multiply bottom, right together
                      call MultiplyMatMat(UmuUnu, Umu_x, Unu_xmu)
                      !# Multiply left, top together, take dagger
                      call MultiplyMatdagMatdag(UmudagUnudag, Umu_xnu, Unu_x)
                      !# multiply two halves together, take trace
                      call RealTraceMultMatMat(P, UmuUnu, UmudagUnudag)
                      sumTrP = sumTrP + P
                      nP = nP + 1
                   end do
                end do
             end do
          end do
       end do
    end do
    call cpu_time(end)
    time = end - start
  end subroutine plaquette


  subroutine plaquette_kolaorder(data, muStart, muEnd, nuEnd, sumTrP, nP, time)
    complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in)  :: data
    integer,                                          intent(in)  :: muStart, muEnd, nuEnd
    real(kind=WP),                                    intent(out) :: sumTrP, time
    integer,                                          intent(out) :: nP
    !"""
    !Calculates the plaquette over muStart to muEnd
    !data is [colour, colour, mu, nx, ny, nz, nt] complex
    !the plaquette over all lattice is muStart=1, muEnd=4, nuEnd=4
    !the spatial plaquette is muStart=1, muEnd=3, nuEnd=3
    !the temporal plaquette is muStart=4, muEnd=5, nuEnd=1
    !returns the sum of plaquettes, number of plaquettes measured,
    !the average plaquette and the time taken to calculate it
    !"""
    integer, dimension(7) :: dataShape
    integer, dimension(4) :: muCoord, nuCoord, coordBase, coord, dataST
    integer :: mu, nu, nx, ny, nz, nt, cc, inc  ! Counters
    ! For intermediate calculating plaquette
    complex(kind=WC), dimension(3, 3) :: Umu_x, Unu_xmu, UmuUnu
    complex(kind=WC), dimension(3, 3) :: Umu_xnu, Unu_x, UmudagUnudag
    real(kind=WP) :: P
    ! Timers
    real :: start, end
    call cpu_time(start)
    dataShape = shape(data)
    dataST = (/ dataShape(4), dataShape(5), dataShape(6), dataShape(7) /)
    !# hold the sum
    sumTrP = 0.0
    !# hold the number measured
    nP = 0
    ! increment nu loop
    inc = 1
    if (muEnd > 4) then
       inc = -1
    end if
    do mu = muStart, muEnd
       if (mu > 4) cycle
       muCoord(:) = 0
       !# This is the shift in mu
       muCoord(mu) = 1
       do nu = mu + 1, nuEnd, inc
          nuCoord(:) = 0
          if (nu .ge. 4 .and. muEnd > 4) cycle
          !# This is the shift in nu
          nuCoord(nu) = 1
          do nx = 1, dataST(1)
             do ny = 1, dataST(2)
                do nz = 1, dataST(3)
                   do nt = 1, dataST(4)
                      !# U_mu(x)
                      coordBase = (/ nx, ny, nz, nt /)
                      coord = coordBase
                      Umu_x = data(:, :, mu, coord(1), coord(2), coord(3), coord(4))
                      !# U_nu(x + amu)
                      coord = coordBase + muCoord
                      !# respect periodic boundary conditions
                      do cc = 1, size(coord)
                         if (coord(cc) > dataST(cc)) then
                            coord(cc) = 1
                         end if
                      end do
                      Unu_xmu = data(:, :, nu, coord(1), coord(2), coord(3), coord(4))
                      !# U_mu(x + anu)
                      coord = coordBase + nuCoord
                      do cc = 1, size(coord)
                         if (coord(cc) > dataST(cc)) then
                            coord(cc) = 1
                         end if
                      end do
                      Umu_xnu = data(:, :, mu, coord(1), coord(2), coord(3), coord(4))
                      !# U_nu(x)
                      coord = coordBase
                      Unu_x = data(:, :, nu, coord(1), coord(2), coord(3), coord(4))
                      !# Multiply bottom, right together
                      call MultiplyMatMat(UmuUnu, Umu_x, Unu_xmu)
                      !# Multiply left, top together, take dagger
                      call MultiplyMatdagMatdag(UmudagUnudag, Umu_xnu, Unu_x)
                      !# multiply two halves together, take trace
                      call RealTraceMultMatMat(P, UmuUnu, UmudagUnudag)
                      sumTrP = sumTrP + P
                      nP = nP + 1
                   end do
                end do
             end do
          end do
       end do
    end do
    call cpu_time(end)
    time = end - start
  end subroutine plaquette_KOLAOrder

  subroutine poly_kolaorder(data, NS, NT,sumTrP, nP, time)
    integer,                                          intent(in)  :: NS, NT
    complex(kind=WC), dimension(3,3, 4, NS, NS, NS, NT), intent(in)  :: data
    real(kind=WP),                                    intent(out) :: sumTrP, time
    integer,                                          intent(out) :: nP
    !"""
    !Calculates the polyakov loop
    !data is [colour, colour, mu, nx, ny, nz, nt] complex
    !returns the sum of polyakov, number of loops measured,
    !the average polyakov and the time taken to calculate it
    !"""
    integer, dimension(7) :: dataShape
    integer, dimension(4) :: coord
    integer :: nx, ny, nz, it  ! Counters
    complex(kind=WC), dimension(3, 3) :: Unu_x, Unu_xt, Unu_Temp
    real(kind=WP) :: P
    ! Timers
    real :: start, end
    call cpu_time(start)
    dataShape = shape(data)
    !# hold the sum
    sumTrP = 0.0
    !# hold the number measured
    nP = 0
    !# loop over all sites
    do nx = 1, NS !dataShape(4)
       do ny = 1, NS ! dataShape(5)
          do nz = 1, NS !dataShape(6)
             ! Get first link U_t(x, 1)
             coord = (/ nx, ny, nz, 1 /)
             Unu_x = data(:, :, 4, coord(1), coord(2), coord(3), coord(4))
             do it = 2, NT - 1 !dataShape(7) - 1
                ! Get middle links U_t(x, t)
                coord = (/ nx, ny, nz, it /)
                Unu_xt = data(:, :, 4, coord(1), coord(2), coord(3), coord(4))
                call MultiplyMatMat(Unu_Temp, Unu_x, Unu_xt)
                Unu_x = Unu_Temp
             end do
             ! get final link U_t(x, NT)
             coord = (/ nx, ny, nz, dataShape(7) /)
             Unu_x = data(:, :, 4, coord(1), coord(2), coord(3), coord(4))
             ! Multiply and trace
             call RealTraceMultMatMat(P, Unu_Temp, Unu_x)
             call MultiplyMatMat(Unu_xt, Unu_Temp, Unu_x)
             sumTrP = sumTrP + P
             nP = nP + 1
          end do
       end do
    end do
    call cpu_time(end)
    time = end - start
  end subroutine poly_kolaorder



  subroutine polyakov(data, sumTrP, nP, time)
    complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in)  :: data
    real(kind=WP),                                    intent(out) :: sumTrP, time
    integer,                                          intent(out) :: nP
    !"""
    !Calculates the polyakov loop
    !data is [nt, nx, ny, nz, mu, colour, colour] complex
    !returns the sum of polyakov, number of loops measured,
    !the average polyakov and the time taken to calculate it
    !"""
    integer, dimension(7) :: dataShape
    integer, dimension(4) :: coord
    integer :: nx, ny, nz, nt  ! Counters
    complex(kind=WC), dimension(3, 3) :: Unu_x, Unu_xt, Unu_Temp
    real(kind=WP) :: P
    ! Timers
    real :: start, end
    call cpu_time(start)
    dataShape = shape(data)
    !# hold the sum
    sumTrP = 0.0
    !# hold the number measured
    nP = 0
    !# loop over all sites
    do nx = 1, dataShape(2)
       do ny = 1, dataShape(3)
          do nz = 1, dataShape(4)
             ! Get first link U_t(x, 1)
             coord = (/ 1, nx, ny, nz /)
             Unu_x = data(coord(1), coord(2), coord(3), coord(4), 1, :, :)
             do nt = 2, dataShape(1) - 1
                ! Get middle links U_t(x, t)
                coord = (/ nt, nx, ny, nz /)
                Unu_xt = data(coord(1), coord(2), coord(3), coord(4), 1, :, :)
                call MultiplyMatMat(Unu_Temp, Unu_x, Unu_xt)
                Unu_x = Unu_Temp
             end do
             ! get final link U_t(x, NT)
             coord = (/ dataShape(1), nx, ny, nz /)
             Unu_x = data(coord(1), coord(2), coord(3), coord(4), 1, :, :)
             ! Multiply and trace
             call RealTraceMultMatMat(P, Unu_Temp, Unu_x)
             call MultiplyMatMat(Unu_xt, Unu_Temp, Unu_x)
             sumTrP = sumTrP + P
             nP = nP + 1
          end do
       end do
    end do
    call cpu_time(end)
    time = end - start
  end subroutine polyakov

end module wloops
