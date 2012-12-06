subroutine stop
USE module_input_parameters,ONLY: MaxLpHaloUsed,MaxMpHaloUsed,mype,lps,lpe,mps,mpe,nprocs,parallelBuild
implicit none
include "gptl.inc"
integer       :: nregions             ! number of gptl regions
integer       :: n,ret                ! index over regions
integer       :: p                    ! index over processors
integer       :: MAXlpHalo            ! Max (over all PEs) lp halo size used
integer       :: MAXmpHalo            ! Max (over all PEs) mp halo size used
integer       :: status
real*8        :: value
real          :: valuemin4, valuemax4 ! For SMS (4-byte quantities)
character(64) :: name                 ! region name
character(80) :: string

! Print timing results to file named timing.mype
ret = gptlpr (mype)

if(parallelBuild) then
  print*
  MAXlpHalo = MaxLpHaloUsed
  MAXmpHalo = MaxMpHaloUsed
!SMS$REDUCE(MAXlpHalo,max)
!SMS$REDUCE(MAXmpHalo,max)
  print"(' Maximum of MAXlpHalo over all processors: ',i0)",MAXlpHalo
  print"(' Maximum of MAXmpHalo over all processors: ',i0)",MAXmpHalo
  print*
  print*,'   PE   lps   lpe   mps   mpe   MAXlpHalo   MAXmpHalo'
!SMS$IGNORE BEGIN
  write(string,100) mype,lps,lpe,mps,mpe,MaxLpHaloUsed,MaxMpHaloUsed
100 format(5i6,i8,i12)
!SMS$IGNORE END
!SMS$INSERT call SMSPrintModeOrdered(string)
endif

print*
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

!ret = gptlprint_memusage ('Memory usage:')

print*,'IPE completed successfully'

end subroutine stop
