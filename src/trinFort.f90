module jclub
  use colours
  use flowAna
  use splineOperations
  use jack
  use openQCDIO
  use wloops
  use types
  implicit none
  !private

  !public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, jclub!"
  end subroutine say_hello
end module jclub
