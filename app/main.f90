program main
   use jclub, only: ReadGaugeFieldUNOpenQCD, polyakov, plaquette, WC, WP
   implicit none
   ! GaugeField params
   !integer, parameter                             :: NS=24, NT=8
   !character(len=128), parameter                  :: gFFile = './Gen2_8x24n7'
   integer, parameter :: NS = 32, NT = 128
   character(len=256), parameter :: gFFile = &
                                    '/home/ryan/Documents/2024/conf/Gen2L/128x32/Gen2l_128x32n1'
   complex(kind=WC), dimension(NT, NS, NS, NS, 4, 3, 3) :: U
   complex(kind=WC), dimension(3, 3, 4, NS, NS, NS, NT) :: UC
   ! plaquette params
   real(kind=WP) :: sumTrP, time
   integer :: nP

   write (*, *) 'NT,NX,NY,NZ,4,3,3'
   U = ReadGaugeFieldUNOpenQCD(gfFile, NS, NS, NS, NT)
   call polyakov(U, sumTrP, np, time)
   write (*, *) 'U Polyakov for', trim(gfFile), ' is ', sumTrp / nP, 'and took', time, 'seconds'
   call plaquette(U, 1, 4, 4, sumTrP, nP, time)
   write (*, *) 'U Plaquette for', trim(gfFile), ' is ', sumTrp / nP, 'and took', time, 'seconds'

   !write(*,*)

   !write(*,*) '3,3,4,NX,NY,NZ,NT'
   !UC = ReadGaugeFieldUNOpenQCDUNKolaOrder(gfFile, NS, NS, NS, NT)
   !call poly_KolaOrder(UC, NS, NT, sumTrP, np, time)
   !write(*,*) 'UC Polyakov for', trim(gfFile), ' is ', sumTrp/nP, 'and took', time, 'seconds'
   !call plaquette_KolaOrder(UC, 1, 4, 4, sumTrP, nP, time)
   !write(*,*) 'UC Plaquette for', trim(gfFile), ' is ', sumTrp/nP, 'and took', time, 'seconds'

end program main
