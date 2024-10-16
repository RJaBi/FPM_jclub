program w0RE
  use jclub, only: WP, complement, setupSplineParams, nx, kx, tx, bcoef, fval, xval, nxv, w1_1d,extrap, ff, splineXDerivMinFunc, setSplineTarget, minFunc, finaliseSplineParams, JackKnife_wp, splineMinFunc
  use stdlib_str2num, only: to_num
  use stdlib_strings, only: find, slice, replace_all, chomp
  use stdlib_math, only: linspace
  use stdlib_io_npy, only: save_npy
  use tomlf, only : toml_table, toml_load, toml_array, get_value, toml_path, len
  use csv_module
  use pyplot_module
  use bspline_module
  use root_module
  implicit none
  ! Input params
  type raggedIntArr
     integer, dimension(:), allocatable :: rag
  end type raggedIntArr
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
  ! csv vars
  character(len=256) :: csvFullPath, csvFileName
  type(csv_file) :: csvFile
  real(kind=WP), dimension(:), allocatable :: flowTime, readCSVArr
  real(kind=WP), dimension(:, :, :), allocatable :: gact4i, gactij!, topcharge, ! flowTime, xi, icon
  logical :: status_ok
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


  ! Load the data
  ! This is awful
  do xx=1, size(xiList)
     icon = 1
     do aa = 1, size(runName)
        xiBase = replace_all(xiPath, "RE", trim(xiList(xx)))
        anaFlow = replace_all(thisFlow, "XI", trim(xiList(xx)))
        anaFlow = replace_all(anaFlow, "NAME", trim(runName(aa)))
        do ii=1, size(iconList(aa)%rag)
           ! Convert into to string
           write(iconStr, '(i0)') iconList(aa)%rag(ii)
           ! write(*,*) trim(iconStr), iconList(aa)%rag(ii)
           csvFileName = replace_all(anaFlow, "ICON", trim(iconStr))
           csvFullPath = trim(xiBase) // '/' // trim(csvFileName)
           ! Read the csv
           call csvFile%read(csvFullPath, header_row=1,status_ok=status_ok, delimiter=' ')
           ! Only read flow-time once
           if (aa == 1 .and. xx == 1 .and. ii == 1) then
              call csvFile%get(2,flowTime,status_ok)
              allocate(gact4i(size(flowTime), size(xiList), ncon), gactij(size(flowTime), size(xiList), ncon))
           end if
           ! Load the data and assign it
           ! write(*,*) csvFullPath
           call csvFile%get(3, readCSVArr, status_ok)
           gact4i(:, xx, icon) = readCSVArr
           deallocate(readCSVArr)
           call csvFile%get(4, readCSVArr, status_ok)
           gactij(:, xx, icon) = readCSVArr
           call csvFile%destroy()
           deallocate(readCSVArr)
           icon = icon + 1
        end do
     end do
  end do

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



contains


  subroutine constructIconList(iconList, iStart, iEnd, iSkip, iSep)
    !! Constructs a list of numbers from iStart to iEnd, excluding
    !! any numbers in iSkip, separtaed by iSep
    integer, dimension(:), allocatable, intent(out) :: iconList
    integer, intent(in) :: iStart, iEnd, iSep
    type(raggedIntArr), intent(in) :: iSkip
    ! counters
    integer :: ii, ncon, icon, jj
    ncon = 0
    iconLoop: do icon = iStart, iEnd, iSep
       do jj =1, size(iSkip%rag)
          if (iSkip%rag(jj) == icon) cycle iconLoop
       end do
       ncon = ncon + 1
    end do iconLoop
    ! Now allocate
    allocate(iconList(ncon))
    ! now assign
    ncon = 1
    iconLoop2: do icon = iStart, iEnd, iSep
       do jj =1, size(iSkip%rag)
          if (iSkip%rag(jj) == icon) cycle iconLoop2
       end do
       iconList(ncon) = icon
       ncon = ncon + 1
    end do iconLoop2
  end subroutine constructIconList


  subroutine processToml(tomlName, producer, eps, tMax, anaDir, xiPath, xiList, xiNumList, runName, iStart, iEnd, iSkip, iSep, targws, targwt, w0PhysMean, w0PhysErr)
    ! Get all the values for the toml
    ! this toml is so much harder than doing it in python...
    character(len=*),  intent(in) :: tomlName
    character(len=:), allocatable, intent(out) :: producer, eps, tMax, anaDir, xiPath
    character(len=128), dimension(:), allocatable, intent(out) :: xiList, runName
    real(kind=WP), dimension(:), allocatable, intent(out) :: xiNumList
    integer, dimension(:), allocatable, intent(out) :: iStart, iEnd, iSep
    type(raggedIntArr), dimension(:), allocatable, intent(out) :: iSkip
    real(kind=WP), intent(out) :: targws, targwt, w0PhysMean, w0PhysErr
    ! working vars
    character(len=:), allocatable :: w0PhysStr
    real(kind=WP) :: w0ErrExp
    ! toml vars
    type(toml_table), allocatable       :: table
    type(toml_array), pointer :: top_array, nested_array
    character(len=:), allocatable :: strRead
    integer :: ii, data_len, jj

    call toml_load(table, tomlName)
    ! Get singleton variables
    call get_value(table, toml_path("data", "xiPath"), xiPath)
    call get_value(table, toml_path("data", "producer"), producer)
    call get_value(table, toml_path("data", "eps"), eps)
    call get_value(table, toml_path("data", "tMax"), tMax)
    call get_value(table, toml_path("analysis", "anaDir"), anaDir)
    call get_value(table, toml_path("analysis", "targws"), targws)
    call get_value(table, toml_path("analysis", "targwt"), targwt)
    ! This first one gets it if it's a string
    call get_value(table, toml_path("analysis", "w0Phys"), w0PhysStr)
    ! This second gets it if its a float
    w0PhysMean = 0.0_WP
    call get_value(table, toml_path("analysis", "w0Phys"), w0PhysMean)
    ! Get array variables
    ! xiList
    call get_value(table, toml_path("data", "xiList"), top_array)
    data_len = len(top_array)
    allocate(xiList(data_len), xiNumList(data_len))
    do ii=1, data_len
       call get_value(top_array, ii, strRead)
       xiList(ii) = strRead
       deallocate(strRead)
    end do
    ! Make an array of the xi in numbers
    xiNumList = to_num(xiList, xiNumList)
    ! runName
    call get_value(table, toml_path("data", "runName"), top_array)
    data_len = len(top_array)
    allocate(runName(data_len), iStart(data_len), iEnd(data_len), iSkip(data_len), iSep(data_len))
    do ii=1, data_len
       call get_value(top_array, ii, strRead)
       runName(ii) = strRead
       deallocate(strRead)
    end do
    ! iStart
    call get_value(table, toml_path("data", "iStart"), top_array)
    do ii=1, len(top_array)
       call get_value(top_array, ii, iStart(ii))
    end do
    ! iEnd
    call get_value(table, toml_path("data", "iEnd"), top_array)
    do ii=1, len(top_array)
       call get_value(top_array, ii, iEnd(ii))
    end do
    ! iSep
    iSep = 1
    call get_value(table, toml_path("data", "iSep"), top_array)
    do ii=1, len(top_array)
       call get_value(top_array, ii, iSep(ii))
    end do
    ! iSkip = 1
    call get_value(table, toml_path("data", "iSkip"), top_array)
    do ii=1, len(top_array)
       call get_value(top_array, ii, nested_array)
       allocate(iSkip(ii)%rag(len(nested_array)))
       iSkip(ii)%rag = -99 ! -99th conf makes no sense of course
       do jj=1, len(nested_array)
          call get_value(nested_array, jj, iSkip(ii)%rag(jj))
       end do
    end do

    if (w0PhysMean == 0.0_WP) then
       if ( find(w0PhysStr, "(") > 0 ) then
          ! Get mean value
          w0PhysMean = to_num(chomp(w0PhysStr, slice(w0PhysStr, find(w0PhysStr, "(") )), w0PhysMean)
          ! Get number for the error (normallised)
          w0PhysErr = to_num(replace_all(replace_all(slice(w0PhysStr, find(w0PhysStr, "(") ), "(", ""), ")", ""), w0PhysErr)
          ! Now need to get 'exponent' of w0PhysErr
          w0ErrExp = -1 + find(w0PhysStr, "(") - find(w0PhysStr, ".")
          w0PhysErr = w0PhysErr * 10**(-w0ErrExp)
       else
          w0PhysErr = 0.0_WP
          w0PhysMean = to_num(w0PhysStr, w0PhysMean)
       end if
    end if

end subroutine processToml



end program w0RE
