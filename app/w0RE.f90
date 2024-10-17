program w0RE
  !  use jclub, only: WP, complement, setupSplineParams, nx, kx, tx, bcoef, fval, xval, nxv, w1_1d,extrap, ff, splineXDerivMinFunc, setSplineTarget, minFunc, finaliseSplineParams, JackKnife_wp, splineMinFunc, constructIconList, raggedIntArr, processToml, loadData, xigCalcs
  use jclub, only: WP, complement, constructIconList, raggedIntArr, processToml, loadData, w0Calcs, JackKnife_wp, col_red, col_blue, col_black, xigCalcs, col_green, spacingCalcs
  use stdlib_strings, only: replace_all
  use stdlib_io_npy, only: save_npy
  use stdlib_math, only: linspace
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
  real(kind=WP), dimension(:,:,:), allocatable :: wijSplineEval, w4iSplineEval ! t, xi, 0:ncon
  real(kind=WP), dimension(:), allocatable :: plotXVal, plotVal, plotErr ! t
  real(kind=WP), dimension(:,:), allocatable :: RESplineEval, a_sSplineEval
  character(len=128) :: plotStr, plotStr2
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
  call loadData(xiList, xiPath, thisFlow, runName, ncon, iconList, flowTime, gact4i, gactij)
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
  ! Do calculations for w0
  call w0Calcs(flowTime, JE4i, JEij, xiNumList, size(flowTime)*2, targws, plotXVal, wijSplineEval, w4iSplineEval, flowTimeForW0, w0ij, w04i)
  ! Plot the W_{ij/4i} data
  allocate(plotVal(size(plotXVal)), plotErr(size(plotXVal)))
  do xx = 1, size(xiList)
     call plt%initialize(grid=.true.,xlabel='$\\tau / a_s$', &
          legend=.true.)
     ! First do wij
     do ii=1, size(plotXVal)
        call JackKnife_wp(ncon, wijSplineEval(ii, xx, :), plotErr(ii))
     end do
     ! plot mean
     plotVal = wijSplineEval(:, xx, 0)
     call plt%add_plot(plotXVal, plotVal, &
          label = '$W_{ij}$', &
          linestyle='--', markersize=0, linewidth=2, istat=istat, color=col_blue)
     ! plot mean + err
     plotVal = wijSplineEval(:, xx, 0) + plotErr
     call plt%add_plot(plotXVal, plotVal, &
          linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_blue)
     ! plot mean - err
     plotVal = wijSplineEval(:, xx, 0) - plotErr
     call plt%add_plot(plotXVal, plotVal, &
          linestyle=':', markersize=0, linewidth=2, istat=istat, label='',color=col_blue)
     call plt%savefig('w0RETest_'//trim(xiList(xx))//'.pdf',istat=istat, pyfile='w0RETest_'//trim(xiList(xx))//'.py')
     ! Then do w4i
     do ii=1, size(plotXVal)
        call JackKnife_wp(ncon, w4iSplineEval(ii, xx, :), plotErr(ii))
     end do
     ! plot mean
     plotVal = w4iSplineEval(:, xx, 0)
     !write(*,*) w4iSplineEval(:, xx, 1)
     !stop
     call plt%add_plot(plotXVal, plotVal, &
          label = '$W_{4i}$', &
          linestyle='--', markersize=0, linewidth=2, istat=istat, color=col_red)
     ! plot mean + err
     plotVal = w4iSplineEval(:, xx, 0) + plotErr
     call plt%add_plot(plotXVal, plotVal, &
          linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
     ! plot mean - err
     plotVal = w4iSplineEval(:, xx, 0) - plotErr
     call plt%add_plot(plotXVal, plotVal, &
          linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
     ! plot target
     plotVal = targws
     call plt%add_plot(plotXVal, plotVal, &
          linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_black)
     call plt%savefig('w0RETest_'//trim(xiList(xx))//'.pdf',istat=istat, pyfile='w0RETest_'//trim(xiList(xx))//'.py')
  end do


  deallocate(plotXVal)
  allocate(RE(size(xiList), 0:ncon))
  RE = w0ij / w04i
  call xigCalcs(RE, xiNumList, size(flowTime)*2, plotXVal, RESplineEval, xig)
  call JackKnife_wp(ncon, xig, err=val)
  write(*,*) 'xig is ', xig(0), ' +- ', val

  ! Plot
  call plt%initialize(grid=.true.,xlabel='$\\xi_{in}$', &
       legend=.true.)
  do ii=1, size(plotXVal)
     call JackKnife_wp(ncon, RESplineEval(ii, :), plotErr(ii))
  end do
  ! plot the spline
  plotVal = RESplineEval(:, 0)
  call plt%add_plot(plotXVal, plotVal, &
       linestyle='-', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
  ! plot mean + err
  plotVal = RESplineEval(:, 0) + plotErr
  call plt%add_plot(plotXVal, plotVal, &
       linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
  ! plot mean - err
  plotVal = RESplineEval(:, 0) - plotErr
  call plt%add_plot(plotXVal, plotVal, &
       linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
  ! Now do the points
  deallocate(plotErr)
  allocate(plotErr(size(xiNumList)))
  do xx=1, size(xiNumList)
     call Jackknife_wp(ncon, RE(xx, :), plotErr(xx))
  end do
  call plt%add_errorbar(xiNumList, RE(:, 0), label='', color=col_blue, istat=istat, linestyle='o',markersize=4,linewidth=0, yerr=plotErr)
  plotVal = 1.0_WP
  call plt%add_plot(plotXVal, plotVal, &
       linestyle='--', markersize=0, linewidth=2, istat=istat, label='$target$', color=col_black)
  plotXVal = xig(0)
  plotVal = linspace(minval(RE(:,0)), maxval(RE(:,0)), size(plotXVal))
  write(plotStr, '(f10.6,a,f7.6)') xig(0), '$ +- $0', val
  call plt%add_plot(plotXVal, plotVal, &
       linestyle='-', markersize=0, linewidth=2, istat=istat, label='$\\xi_{g} = '//trim(plotStr)//'$', color=col_green)
    call plt%add_plot(plotXVal + val, plotVal, &
         linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_green)
    call plt%add_plot(plotXVal - val, plotVal, &
       linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_green)
  call plt%savefig('w0RETest_RE.pdf',istat=istat, pyfile='w0RETest_RE.py')
  ! Do the lattice spacing
  deallocate(plotXVal)
  call spacingCalcs(flowTimeForW0, xiNumList, xig, size(flowTime)*2, w0PhysMean, w0PhysErr, a_sSplineEval, plotXVal, a_s, a_sSys)

  call Jackknife_wp(ncon, a_s, val)
  call plt%initialize(grid=.true.,xlabel='$\\xi_{in}$', &
       legend=.true.)



  ! Plot the spline
  deallocate(plotErr)
  allocate(plotErr(size(plotXVal)))
    do ii=1, size(plotXVal)
     call JackKnife_wp(ncon, a_sSplineEval(ii, :), plotErr(ii))
  end do
  ! plot the spline
  plotVal = a_sSplineEval(:, 0)
  call plt%add_plot(plotXVal, plotVal, &
       linestyle='-', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
  ! plot mean + err
  plotVal = a_sSplineEval(:, 0) + plotErr
  call plt%add_plot(plotXVal, plotVal, &
       linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
  ! plot mean - err
  plotVal = a_sSplineEval(:, 0) - plotErr
  call plt%add_plot(plotXVal, plotVal, &
       linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
  ! Plot the lattice spacing horizontally
  plotXVal = linspace(minval(xiNumList), maxval(xiNumList), size(plotXVal))
  plotVal = a_s(0)
  write(plotStr, '(f10.6,a,f7.6,a,f7.6,a,f7.6,a)') a_s(0), '(0', val, ')(0',a_sSys,')[0',(val**2.0_WP + a_sSys**2.0_WP)**0.5_WP,']'
  val = (val**2.0_WP + a_sSys**2.0_WP)**0.5_WP
  call plt%add_plot(plotXVal, plotVal, &
       linestyle='-', markersize=0, linewidth=2, istat=istat, label='$a_{s} = '//trim(plotStr)//'$', color=col_green)
  call plt%add_plot(plotXVal + val, plotVal, &
       linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_green)
  call plt%add_plot(plotXVal - val, plotVal, &
       linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_green)

  ! Plot the anisotropy vertically
  plotXVal = xig(0)
  !plotVal = linspace(minval(RE(:,0)), maxval(RE(:,0)), size(plotXVal))
  plotVal = linspace(minval(a_sSplineEval(:,0)), maxval(a_sSplineEval(:,0)), size(plotXVal))
  write(plotStr, '(f10.6,a,f7.6)') xig(0), '$ +- $0', val
  call plt%add_plot(plotXVal, plotVal, &
       linestyle='-', markersize=0, linewidth=2, istat=istat, label='$\\xi_{g} = '//trim(plotStr)//'$', color=col_black)
    call plt%add_plot(plotXVal + val, plotVal, &
         linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_black)
    call plt%add_plot(plotXVal - val, plotVal, &
       linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_black)

    call plt%savefig('w0RETest_as.pdf',istat=istat, pyfile='w0RETest_as.py')
  ! Now do the points

  ! write(*,*) a_s(0), val, a_sSys, (val**2.0_WP + a_sSys**2.0_WP)**0.5_WP


  deallocate(runName, xiList, xiNumList)
  deallocate(iStart, iEnd, iSep, iSkip)
  deallocate(iconList)
  deallocate(flowTime, gact4i, gactij)
  deallocate(JE4i, JEij)
  ! more?

  write(*,*) 'Done'


end program w0RE
