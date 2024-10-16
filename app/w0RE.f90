program w0RE
  use jclub, only: WP, complement, setupSplineParams, nx, kx, tx, bcoef, fval, xval, nxv, w1_1d,extrap, ff, splineXDerivMinFunc, setSplineTarget, minFunc, finaliseSplineParams, JackKnife_wp, splineMinFunc, constructIconList, raggedIntArr, processToml, loadData
  use stdlib_strings, only: replace_all
  use stdlib_math, only: linspace
  use stdlib_io_npy, only: save_npy
  use pyplot_module
  use bspline_module
  use root_module
  implicit none
  ! Input params
  character(len=:), allocatable :: producer, eps, tMax, anaDir, w0PhysStr, xiPath
  character(len=128), dimension(:), allocatable :: xiList, runName
  real(kind=WP), dimension(:), allocatable :: xiNumList
  integer, dimension(:), allocatable :: iStart, iEnd, iSep
  type(raggedIntArr), dimension(:), allocatable :: iSkip
  real(kind=WP) :: targws, targwt, w0PhysMean, w0PhysErr
  ! IO vars
  type(raggedIntArr), dimension(:), allocatable :: iconList
  character(len=128) :: xiBase, thisFlow, anaFlow, iconStr
  character(len=128), parameter :: flowBase = 'flow.NAME_xiXI_eEPS_tTMAX.ICON'
  ! counters
  integer :: xx, icon, ncon, aa, ii
  real(kind=WP), dimension(:), allocatable :: flowTime
  real(kind=WP), dimension(:, :, :), allocatable :: gact4i, gactij!, topcharge, ! flowTime, xi, icon
  ! Plotting vars
  type(pyplot) :: plt
  integer :: istat
  ! Jack-knifed vars
  real(kind=WP), dimension(:,:,:), allocatable :: JE4i, JEij ! flowTime, xi, 0:ncon
  real(kind=WP), dimension(:,:), allocatable :: w0ij, w04i, RE! xi, 0:ncon
  real(kind=WP), dimension(:,:), allocatable :: flowTimeForW0 ! xi, 0:ncon
  real(kind=WP), dimension(:), allocatable :: a_s, xig  ! xi, 0:ncon
  ! spline vars
  integer :: iflag, inbvx, iloy
  real(kind=WP) :: val
  ! root finder
  real(kind=WP) :: xr, fr
  integer :: rflag
  ! a_s variables
  real(kind=WP) :: a_sMean, a_sStat, a_sSys
  ! Setup and Get all the parameters from the input toml file
  call processToml('/home/ryan/Documents/2024/Gen2/G2_wflow.toml', producer, eps, tMax, anaDir, xiPath, xiList, xiNumList, runName, iStart, iEnd, iSkip, iSep, targws, targwt, w0PhysMean, w0PhysErr)

  thisFlow = replace_all(flowBase, 'EPS', trim(eps))
  thisFlow = replace_all(thisFlow, 'TMAX', trim(tmax))
  allocate(iconList(size(runName)))
  ncon = 0
  do aa=1, size(runName)
     call constructIconList(iconList(aa)%rag, iStart(aa), iEnd(aa), iSkip(aa), iSep(aa))
     ncon = ncon + size(iconList(aa)%rag)
  end do



  call loadData(xiList, xiPath, thisFlow, runName, ncon, iconList, flowTime, gact4i, gactij)
!csv!
!csv!  ! Load the data
!csv!  ! This is awful
!csv!  do xx=1, size(xiList)
!csv!     icon = 1
!csv!     do aa = 1, size(runName)
!csv!        xiBase = replace_all(xiPath, "RE", trim(xiList(xx)))
!csv!        anaFlow = replace_all(thisFlow, "XI", trim(xiList(xx)))
!csv!        anaFlow = replace_all(anaFlow, "NAME", trim(runName(aa)))
!csv!        do ii=1, size(iconList(aa)%rag)
!csv!           ! Convert into to string
!csv!           write(iconStr, '(i0)') iconList(aa)%rag(ii)
!csv!           ! write(*,*) trim(iconStr), iconList(aa)%rag(ii)
!csv!           csvFileName = replace_all(anaFlow, "ICON", trim(iconStr))
!csv!           csvFullPath = trim(xiBase) // '/' // trim(csvFileName)
!csv!           ! Read the csv
!csv!           call csvFile%read(csvFullPath, header_row=1,status_ok=status_ok, delimiter=' ')
!csv!           ! Only read flow-time once
!csv!           if (aa == 1 .and. xx == 1 .and. ii == 1) then
!csv!              call csvFile%get(2,flowTime,status_ok)
!csv!              allocate(gact4i(size(flowTime), size(xiList), ncon), gactij(size(flowTime), size(xiList), ncon))
!csv!           end if
!csv!           ! Load the data and assign it
!csv!           ! write(*,*) csvFullPath
!csv!           call csvFile%get(3, readCSVArr, status_ok)
!csv!           gact4i(:, xx, icon) = readCSVArr
!csv!           deallocate(readCSVArr)
!csv!           call csvFile%get(4, readCSVArr, status_ok)
!csv!           gactij(:, xx, icon) = readCSVArr
!csv!           call csvFile%destroy()
!csv!           deallocate(readCSVArr)
!csv!           icon = icon + 1
!csv!        end do
!csv!     end do
!csv!  end do
!csv!
  ! Do Jackknifes
  allocate(JE4i(size(flowTime), size(xiList), 0:ncon), JEij(size(flowTime), size(xiList), 0:ncon))
  do ii=1, size(flowTime)
     do xx=1, size(xiList)
        ! Take Mean
        call complement(JE4i(ii, xx, 0), gact4i(ii, xx, :) * flowTime(ii)**2.0)
        call complement(JEij(ii, xx, 0), gactij(ii, xx, :) * flowTime(ii)**2.0)
        ! Take 1st order jackknifes
        call complement(ncon, JE4i(ii, xx, 1:), gact4i(ii, xx, :)* flowTime(ii)**2.0)
        call complement(ncon, JEij(ii, xx, 1:), gactij(ii, xx, :)* flowTime(ii)**2.0)
     end do
  end do



  call plt%initialize(grid=.true.,xlabel='$\\tau / a_s$', &
       legend=.true.)

  call setupSplineParams(size(flowTime), 4, size(flowTime)*2)
  allocate(w0ij(size(xiList), 0:ncon), w04i(size(xiList), 0:ncon), RE(size(xiList), 0:ncon))
  allocate(flowTimeForW0(size(xiList), 0:ncon))
  do xx=1, size(xiList)
     do icon=0, ncon
        !Fit a spline to calculate w0
        !have to set these before the first evaluate call:
        inbvx = 1
        iloy  = 1
        ! initialize and fit spline
        ! for wij
        call db1ink(flowTime,nx,JEij(:, xx, icon),kx,0,tx,bcoef,iflag)
        call setSplineTarget(targws)
        ff=>splineXDerivMinFunc
        ! Match python tolerance better
        call root_scalar('brent', minFunc, minval(flowTime), maxval(flowTime), xr, fr, rflag, rtol=0.0000000000000001_WP)
        !write(*,*) 'fr(',xr,') = ', fr
        flowTimeForW0(xx, icon) = xr**0.5_WP
        w0ij(xx, icon) = fr + targws
        ! for w4i
        !have to set these before the first evaluate call:
        inbvx = 1
        iloy  = 1
        call db1ink(flowTime,nx,JE4i(:, xx, icon),kx,0,tx,bcoef,iflag)
        call db1val(flowTimeForW0(xx, icon) **2.0_WP,1,tx,nx,kx,bcoef,w04i(xx, icon),iflag,inbvx,w1_1d,extrap=extrap)
        w04i(xx, icon) = w04i(xx, icon) * (flowTimeForW0(xx, icon) **2.0_WP) * (xiNumList(xx)**2.0)
        ! write(*,*) trim(xiList(xx)), w0ij(xx, icon), w04i(xx, icon), w0ij(xx, icon) / w04i(xx, icon)
        RE(xx, icon) = w0ij(xx, icon) / w04i(xx, icon)
     end do
  end do


  ! Now we fit a spline to the RE data
  !do xx = 1, size(xiList)
     !call JackKnife_wp(ncon, RE(xx, :), val)
     !write(*,*) trim(xiList(xx)), RE(xx, 0), val
     !call JackKnife_wp(ncon, flowTimeForW0(xx, :), val)
     !write(*,*) trim(xiList(xx)), flowTimeForW0(xx, 0), val
     !call JackKnife_wp(ncon, w0ij(xx, :) - 0.15_WP, val)
     !write(*,*) trim(xiList(xx)), w0ij(xx, 0) - 0.15_WP, val
     !call JackKnife_wp(ncon, w04i(xx, :), val)
     !write(*,*) trim(xiList(xx)), w04i(xx, 0), val
  !end do

  allocate(xig(0:ncon))

  ! finalise the old splines
  call finaliseSplineParams
  ! setup the new ones
  call setupSplineParams(size(xiList), 4, size(flowTime)*2)
  do icon=0, ncon
     inbvx = 1
     iloy  = 1
     call db1ink(xiNumList,nx,RE(:, icon),kx,0,tx,bcoef,iflag)
     call setSplineTarget(1.0_WP)
     ff=>splineMinFunc
     call root_scalar('brent', minFunc, minval(xiNumList), maxval(xiNumList), xr, fr, rflag, rtol=0.0000000000000001_WP)
     xig(icon) = xr
  end do
  call JackKnife_wp(ncon, xig, err=val)
  write(*,*) 'xig is ', xig(0), ' +- ', val

  call save_npy('xig.npy', xig)

  ! Now for the lattice spacing
  allocate(a_s(0:ncon))
  do icon = 0, ncon
     ! Need to fit flowTimeForW0 as function of xi
     ! initialize and fit spline
     inbvx = 1
     iloy  = 1
     call db1ink(xiNumList,nx,flowTimeForW0(:, icon),kx,0,tx,bcoef,iflag)
     ! Evaluate spline at xig
     call db1val(xig(icon), 0, tx, nx, kx, bcoef, a_s(icon), iflag, inbvx, w1_1d, extrap=extrap)
  end do

  write(*,*) a_s(0)

  a_s = w0PhysMean / a_s
  a_sMean = a_s(0)
  call JackKnife_wp(ncon, a_s, err=a_sStat)
  !write(*,*) 'a_s (Stat) is ', a_s(0), ' +- ', val
  a_s = (w0PhysMean + w0PhysErr)/ (W0PhysMean / a_s)
  a_sSys = abs(a_s(0) - a_sMean)
  !call JackKnife_wp(ncon, a_s, err=val)
  !write(*,*) 'a_s (+) is ', a_s(0), ' +- ', val
  a_s = (w0PhysMean - w0PhysErr)/ ((W0PhysMean + w0PhysErr)/ a_s)
  if (abs(a_s(0) - a_sMean) > a_sSys) then
     a_sSys = abs(a_s(0) - a_sMean)
  end if
  val = (a_sStat**2.0_WP + a_sSys**2.0_WP)**0.5_WP
  write(*,*) 'a_s=', a_sMean, '+- Stat', a_sStat, '+- Sys', a_sSys
  write(*,*) 'a_s=', a_sMean, '+- combined ', val
  !call JackKnife_wp(ncon, a_s, err=val)
  !write(*,*) 'a_s (-) is ', a_s(0), ' +- ', val

  !a_s(:, :) = w0PhysMean

  !  xval = linspace(0.0_WP, maxval(flowTime), nxv)
!  do ii=1,nxv
!     call db1val(xval(ii),1,tx,nx,kx,bcoef,val,iflag,inbvx,w1_1d,extrap=extrap)
!     fval(ii) = (xiNumList(xx)**2.0) * xval(ii) * val  ! save it for plot
!  end do
!  call plt%add_plot(xval, fval, &
!       label='$W_{ij}$',&
!       linestyle='ko',markersize=5,linewidth=2,istat=istat)


  !call plt%savefig('w0RETest.pdf',istat=istat, pyfile='w0ReTest.py')
  !xval = linspace(0.0_WP, maxval(flowTime), nxv)
  !do ii=1,nxv
  !   call db1val(xval(ii),1,tx,nx,kx,bcoef,val,iflag,inbvx,w1_1d,extrap=extrap)
  !   fval(ii) = xval(ii) * val  ! save it for plot
  !end do
  ! Set the target for the spline


  ! Simple plot


  call finaliseSplineParams
  deallocate(runName, xiList, xiNumList)
  deallocate(iStart, iEnd, iSep, iSkip)
  deallocate(iconList)
  deallocate(flowTime, gact4i, gactij)
  deallocate(JE4i, JEij)

  write(*,*) 'Done'


end program w0RE
