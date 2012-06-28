      subroutine stop
      implicit none
      include "gptl.inc"
      integer       :: nregions             ! number of gptl regions
      integer       :: n,ret                ! index over regions
      real*8        :: value
      real          :: valuemin4, valuemax4 ! For SMS (4-byte quantities)
      character(64) :: name                 ! region name

      ret = gptlget_nregions (0, nregions)
      do n=0,nregions-1
        name = ' '
        ret = gptlget_regionname (0, n, name)
        ret = gptlget_wallclock (trim (name), 0, value)
        valuemin4 = value
        valuemax4 = value
!SMS$REDUCE(valuemin4,min)
!SMS$REDUCE(valuemax4,max)
        print"(1x,A20,F13.3,F15.3)", trim (name), valuemin4, valuemax4
      end do

! Print timing results to file named timing.0
      ret = gptlpr (0)

      end subroutine stop
