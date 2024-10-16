module flowAna
  use types, only: WP, raggedIntArr
  use tomlf, only : toml_table, toml_load, toml_array, get_value, toml_path, len
  use stdlib_str2num, only: to_num
  use stdlib_strings, only: find, slice, replace_all, chomp
  use csv_module
  implicit none

  private

  public :: constructIconList, processToml, loadData

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

  subroutine xigCalcs(flowTime, gact4i, gactij, anaDir, xiNumList, wijSplineEval, w4iSplineEval)
    real(kind=WP), dimension(:), intent(in) :: flowTime
    real(kind=WP), dimension(:,:,:), intent(in) :: gact4i, gactij
    character(len=*), intent(in) :: anaDir
    real(kind=WP), dimension(:), intent(in) :: xiNumList
    real(kind=WP), dimension(:,:,:), intent(out) :: wijSplineEval, w4iSplineEval




  end subroutine xigCalcs


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
