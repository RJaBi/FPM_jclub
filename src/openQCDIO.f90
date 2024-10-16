module openqcdio
  use types, only: wc, wp
  implicit none
  !private
  public
  !public :: readgaugefieldUNopenqcd, readgaugefieldUNopenqcdUNkolaorder

  !public :: determinant

contains

  !! maps an integer a to the set of integers [1,b] i.e. positive integers with cycle length b.
  elemental function modc(a,b) result (c)
    !! taken directly from kola
    integer, intent(in) :: a,b
    integer :: c
    !c = a - ((a-1)/b)*b
    c = modulo(a-1,b)+1
  end function modc


  function determinant(matrix) result(det)
    use, intrinsic ::  iso_c_binding, only: c_double_complex, c_double
    complex(c_double_complex), dimension(3,3), intent(in) :: matrix
    complex(c_double_complex) :: det
!
    det = matrix(1,1)*(matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) - &
         matrix(1,2)*(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) + &
         matrix(1,3)*(matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))
  end function determinant

  function readgaugefieldUNopenqcd(filename,nx, ny, nz, nt) result(u_xd)
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: nx, ny, nz, nt
    complex(kind=wc), dimension(nt,nx,ny,nz,4,3,3) :: u_xd
    complex(kind=wc), dimension(:,:,:,:,:,:,:), allocatable :: u_tmp

    integer, parameter :: infl = 107
    ! header info
    real(kind=wp) :: plaq
    integer :: ntdim, nxdim, nydim, nzdim
    ! counters
    integer :: it, ix, iy, iz, mu, id
    integer :: jx, jy, jz, jt
    integer, dimension(4) :: dmu
    open(infl, file=trim(filename), form='unformatted', access='stream', status='old', action='read', convert='little_endian')
    read(infl) ntdim, nxdim, nydim, nzdim, plaq
    write(*,*) ntdim, nxdim, nydim, nzdim, plaq
    ! need the tmp array to make the cshift later work
    allocate(u_tmp(ntdim,nydim,nzdim,nzdim,4,3,3))
    ! z varies quickest, then y, then x, then t
    do it = 1,ntdim; do ix = 1,nxdim; do iy = 1,nydim; do iz = 1,nzdim
       if ( modulo(ix+iy+iz+it-4,2) == 0 ) cycle  ! format only considers odd points

       do id = 1,4
          mu = modc(id-1,4)  ! time dimension first: mu = 4, 1, 2, 3

          dmu(:) = 0
          dmu(mu) = 1

          ! get the backward site under periodic boundary conditions
          jx = modc(ix-dmu(1),nxdim)
          jy = modc(iy-dmu(2),nydim)
          jz = modc(iz-dmu(3),nzdim)
          jt = modc(it-dmu(4),ntdim)
          ! read the forward and backward links in mu direction
          read(infl) u_tmp(it,ix,iy,iz,mu,:,:)
          read(infl) u_tmp(jt,jx,jy,jz,mu,:,:)
          u_tmp(it,ix,iy,iz,mu,:,:) = transpose(u_tmp(it,ix,iy,iz,mu,:,:))
          u_tmp(jt,jx,jy,jz,mu,:,:) = transpose(u_tmp(jt,jx,jy,jz,mu,:,:))

          !call fixsu3matrix(u_g(mu,ix,iy,iz,it))
          !call fixsu3matrix(u_g(mu,jx,jy,jz,jt))
       end do
    end do; end do; end do; end do

    close(infl)

    u_xd = cshift(u_tmp, -1, dim=5)
    deallocate(u_tmp)
  end function readgaugefieldUNopenqcd

    function readgaugefieldUNopenqcdUNkolaorder(filename,nx, ny, nz, nt) result(u_xd)
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: nx, ny, nz, nt
    complex(kind=wc), dimension(3,3,4,nx,ny,nz,nt) :: u_xd
    complex(kind=wc), dimension(:,:,:,:,:,:,:), allocatable :: u_tmp
    complex(kind=wc), dimension(:,:,:,:,:,:,:), allocatable :: u_tc

    integer, parameter :: infl = 107
    ! header info
    real(kind=wp) :: plaq
    integer :: ntdim, nxdim, nydim, nzdim
    ! counters
    integer :: it, ix, iy, iz, mu, id
    integer :: jx, jy, jz, jt
    integer, dimension(4) :: dmu
    open(infl, file=trim(filename), form='unformatted', access='stream', status='old', action='read', convert='little_endian')
    read(infl) ntdim, nxdim, nydim, nzdim, plaq
    ! need the tmp array to make the cshift later work
    allocate(u_tmp(ntdim,nydim,nzdim,nzdim,4,3,3))
    allocate(u_tc(ntdim,nydim,nzdim,nzdim,4,3,3))
    write(*,*) ntdim, nxdim, nydim, nzdim, plaq
    ! z varies quickest, then y, then x, then t
    do it = 1,ntdim; do ix = 1,nxdim; do iy = 1,nydim; do iz = 1,nzdim
       if ( modulo(ix+iy+iz+it-4,2) == 0 ) cycle  ! format only considers odd points

       do id = 1,4
          mu = modc(id-1,4)  ! time dimension first: mu = 4, 1, 2, 3

          dmu(:) = 0
          dmu(mu) = 1

          ! get the backward site under periodic boundary conditions
          jx = modc(ix-dmu(1),nxdim)
          jy = modc(iy-dmu(2),nydim)
          jz = modc(iz-dmu(3),nzdim)
          jt = modc(it-dmu(4),ntdim)
          ! read the forward and backward links in mu direction
          read(infl) u_tmp(it,ix,iy,iz,mu,:,:)
          read(infl) u_tmp(jt,jx,jy,jz,mu,:,:)
          u_tmp(it,ix,iy,iz,mu,:,:) = transpose(u_tmp(it,ix,iy,iz,mu,:,:))
          u_tmp(jt,jx,jy,jz,mu,:,:) = transpose(u_tmp(jt,jx,jy,jz,mu,:,:))

          !call fixsu3matrix(u_g(mu,ix,iy,iz,it))
          !call fixsu3matrix(u_g(mu,jx,jy,jz,jt))
       end do
    end do; end do; end do; end do
    close(infl)
    u_tc = cshift(u_tmp, -1, dim=5)
    ! reshape data manually
    do ix=1, nx
       do iy=1, ny
          do iz=1, nz
             do it=1, nt
                u_xd(:,:,1,ix,iy,iz,it) = (u_tc(it,ix,iy,iz,2,:,:)) ! x
                u_xd(:,:,2,ix,iy,iz,it) = (u_tc(it,ix,iy,iz,3,:,:)) ! y
                u_xd(:,:,3,ix,iy,iz,it) = (u_tc(it,ix,iy,iz,4,:,:)) ! z
                u_xd(:,:,4,ix,iy,iz,it) = (u_tc(it,ix,iy,iz,1,:,:)) ! t
             end do
          end do
       end do
    end do
    deallocate(u_tc, u_tmp)
  end function readgaugefieldUNopenqcdUNkolaorder


end module openqcdio
