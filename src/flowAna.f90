module flowAna
  use types, only: WP, raggedIntArr
  use splineOperations, only: setupSplineParams, nx, kx, tx, bcoef, fval, xval, &
       nxv, w1_1d, extrap, ff, splineXDerivMinFunc, setSplineTarget, minFunc, &
       finaliseSplineParams, splineMinFunc
  use tomlf, only : toml_table, toml_load, toml_array, get_value, toml_path, len
  use stdlib_str2num, only: to_num
  use stdlib_strings, only: find, slice, replace_all, chomp
  use stdlib_math, only: linspace
  use csv_module
  use bspline_module
  use root_module
  use jack, only: jackknife_wp
  implicit none

  private

  public :: constructIconList, processToml, loadData, w0Calcs, xigCalcs, spacingCalcs

contains

  subroutine loadData(xiList, xiPath, thisFlow, runName, ncon, iconList, flowTime, gact4i, gactij)
    ! Loads the data from all the csv's
    character(len=128), dimension(:), intent(in) :: xiList, runName
    character(len=*), intent(in) :: xiPath, thisFlow
    integer, intent(in) :: ncon
    type(raggedIntArr), dimension(:), intent(in) :: iconList
    real(kind=WP), dimension(:), allocatable, intent(out) :: flowTime
    real(kind=WP), dimension(:,:,:), allocatable, intent(out) :: gact4i, gactij
    ! counters
    integer :: icon, xx, aa, ii
    ! strings
    character(len=128) :: xiBase, anaFlow, iconStr, csvFileName, csvFullPath
    ! csv
    type(csv_file) :: csvFile
    logical :: status_ok
    real(kind=WP), dimension(:), allocatable :: readCSVArr
    ! Load the data
    ! This is awful but works
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
  end subroutine loadData


  subroutine spacingCalcs(flowTimeForW0, xiNumList, xig, nEval, w0PhysMean, w0PhysErr, a_sSplineEval, xEval, a_s, a_sSys)
    real(kind=WP), dimension(:, 0:), intent(in) :: flowTimeForW0 ! xi, 0:ncon
    real(kind=WP), dimension(:), intent(in) :: xiNumList ! xi
    real(kind=WP), dimension(0:), intent(in) :: xig ! 0:ncon
    integer, intent(in) :: nEval
    real(kind=WP), intent(in) :: w0PhysMean, w0PhysErr  ! in fm
    real(kind=WP), dimension(:, 0:), intent(out) :: a_sSplineEval  ! nEval, 0:ncon
    real(kind=WP), dimension(:), intent(out) :: xEval ! nEval
    real(kind=WP), dimension(0:), intent(out) :: a_s
    real(kind=WP), intent(out) :: a_sSys
    ! counters mostly
    integer :: icon, ncon, inbvx, iloy, iflag, rflag, ii
    ! evaluated single numbers
    real(kind=WP) :: val
    ncon = size(flowTimeForW0(1,:)) - 1
    call setupSplineParams(size(xiNumList), 4, nEval)
    xEVal = linspace(minval(xiNumList), maxval(xiNumList), nEval)
    do icon = 0, ncon
       ! Need to fit flowTimeForW0 as function of xi
       ! initialize and fit spline
       inbvx = 1
       iloy  = 1
       call db1ink(xiNumList,nx,flowTimeForW0(:, icon),kx,0,tx,bcoef,iflag)
       ! Evaluate spline at xig
       call db1val(xig(icon), 0, tx, nx, kx, bcoef, a_s(icon), iflag, inbvx, w1_1d, extrap=extrap)
       ! do the spline eval
       do ii=1, nEval
          call db1val(xEval(ii),0,tx,nx,kx,bcoef,val,iflag,inbvx,w1_1d,extrap=extrap)
          a_sSplineEval(ii, icon) = W0PhysMean / val
       end do
    end do
    a_s = w0PhysMean / a_s
    val = a_s(0)
    ! Get systematic from uncertainty in w0Phys
    ! Check +
    val = (w0PhysMean + w0PhysErr)/ (W0PhysMean / a_s(0))
    a_sSys = abs(a_s(0) - val)
    ! Check -
    val = (w0PhysMean - w0PhysErr)/ ((W0PhysMean)/ a_s(0))
    if (abs(a_s(0) - val) > a_sSys) then
       a_sSys = abs(a_s(0) - val)
    end if
  call finaliseSplineParams()

  end subroutine spacingCalcs


  subroutine xigCalcs(RE, xiNumList, nEval, xEval, RESplineEval, xig)
    real(kind=WP), dimension(:, 0:), intent(in) :: RE  ! xi, 0:ncon
    real(kind=WP), dimension(:), intent(in) :: xiNumList
    integer, intent(in) :: neVal
    real(kind=WP), dimension(:), intent(out) :: xEval
    real(kind=WP), dimension(:, 0:), intent(out) :: ReSplineEval ! nEval, 0:ncon
    real(kind=WP), dimension(0:), intent(out) :: xig ! 0:ncon
    ! counters mostly
    integer :: icon, ncon, inbvx, iloy, iflag, rflag, ii
    ! evaluated single numbers
    real(kind=WP) :: val, xr, fr
    ncon = size(RE(1,:)) - 1
    ! allocate(xEval(nEval), RESplineEval(nEval,0:ncon))
    xEVal = linspace(minval(xiNumList), maxval(xiNumList), nEval)
    ! allocate(xig(0:ncon))
    call setupSplineParams(size(xiNumList), 4, nEval)
    call setSplineTarget(1.0_WP)
    do icon=0, ncon
       inbvx = 1
       iloy = 1
       call db1ink(xiNumList, nx, RE(:, icon), kx, 0, tx, bcoef, iflag)
       ! get the xig
       ff=>splineMinFunc
       call root_scalar('brent', minFunc, minval(xiNumList), maxval(xiNumList), xr, fr, rflag, rtol=0.0000000000000001_WP)
       xig(icon) = xr
       ! do the spline eval
       do ii=1, nEval
          call db1val(xEval(ii),0,tx,nx,kx,bcoef,val,iflag,inbvx,w1_1d,extrap=extrap)
          RESplineEval(ii, icon) = val
       end do
    end do
    call finaliseSplineParams
  end subroutine xigCalcs

  subroutine w0Calcs(flowTime, JE4i, JEij, xiNumList, nEval, targws, xEval, wijSplineEval, w4iSplineEval, flowTimeForW0, w0ij, w04i)
    real(kind=WP), dimension(:), intent(in) :: flowTime
    real(kind=WP), dimension(:,:,0:), intent(in) :: JE4i, JEij  ! t, xi, 0:ncon
    real(kind=WP), dimension(:), intent(in) :: xiNumList
    integer, intent(in) :: nEval
    real(kind=WP), intent(in) :: targws
    real(kind=WP), dimension(:,:,0:),  intent(out) :: wijSplineEval, w4iSplineEval ! t, xi, 0:ncon
    real(kind=WP), dimension(:, 0:), intent(out) :: flowTimeForW0 ! xi, 0:ncon
    real(kind=WP), dimension(:,0:), intent(out) :: w0ij, w04i ! xi, 0:ncon
    real(kind=WP), dimension(:), intent(out) :: xEval  ! t
    ! counters
    integer :: xx, ncon, icon, ii
    ! spline vars
    integer :: inbvx, iloy, iflag
    real(kind=WP) :: fEval
    ! root finder
    real(kind=WP) :: xr, fr
    integer :: rflag
    ncon = size(JE4i(1,1,:)) - 1
    call setupSplineParams(size(flowTime), 4, nEval)
    !allocate(w0ij(size(xiNumList), 0:ncon), w04i(size(xiNumList), 0:ncon))
    !allocate(flowTimeForW0(size(xiNumList), 0:ncon))
    !allocate(wijSplineEval(nEval, size(xiNumList), 0:ncon), w4iSplineEval(nEval, size(xiNumList), 0:ncon))
    !allocate(xEval(nEval))

    xEval = linspace(0.0_WP, maxval(flowTime), nEval)

    do xx=1, size(xiNumList)
       do icon=0, ncon
          !Fit a spline to calculate w0
          !have to set these before the first evaluate call:
          inbvx = 1
          iloy  = 1
          ! initialize and fit spline to <E(t)*t^2>
          ! for wij
          call db1ink(flowTime,nx,JEij(:, xx, icon),kx,0,tx,bcoef,iflag)
          call setSplineTarget(targws)
          ! Calculate <w(t) = t * d <E(t)*t^2> / dt> and find intercept with targws
          ff=>splineXDerivMinFunc
          call root_scalar('brent', minFunc, minval(flowTime), maxval(flowTime), xr, fr, rflag, rtol=0.0000000000000001_WP)
          flowTimeforW0(xx, icon) = xr**0.5
          w0ij(xx, icon) = fr + targws
          ! Evaluate <w(t) = t * d <E(t)*t^2> / dt> for plot
          do ii=1, nEval
             wijSplineEval(ii, xx, icon) = ff(xEval(ii)) + targws
          end do
          ! Now do temporal
          ! for w4i
          !have to set these before the first evaluate call:
          inbvx = 1
          iloy  = 1
          call db1ink(flowTime,nx,JE4i(:, xx, icon),kx,0,tx,bcoef,iflag)
          call db1val(flowTimeForW0(xx, icon)**2.0_WP,1,tx,nx,kx,bcoef,w04i(xx, icon),iflag,inbvx,w1_1d,extrap=extrap)
          w04i(xx, icon) = w04i(xx, icon) * (flowTimeForW0(xx, icon) **2.0_WP) * (xiNumList(xx)**2.0)
          ! Evaluate xi**2.0 * <w(t) = xi**2.0 * t * d <E(t)*t^2> / dt>for plot
          do ii=1, nEval
             call db1val(xEval(ii),1,tx,nx,kx,bcoef,fEval,iflag,inbvx,w1_1d,extrap=extrap)
             w4iSplineEval(ii, xx, icon) = (xiNumList(xx)**2.0) * xEval(ii) * fEval
          end do
       end do
    end do
    call finaliseSplineParams

  end subroutine w0Calcs


  subroutine processToml(tomlName, producer, eps, tMax, anaDir, xiPath, xiList, &
       xiNumList, runName, iStart, iEnd, iSkip, iSep, targws, targwt, w0PhysMean, w0PhysErr)
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




end module flowAna
