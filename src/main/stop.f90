subroutine stop
USE module_input_parameters,ONLY: MaxLpHaloUsed,MaxMpHaloUsed,mype,lps,lpe,mps,mpe,nprocs,parallelBuild
implicit none
include "gptl.inc"
integer                   :: nregions ! number of gptl regions
integer                   :: n,ret    ! index over regions
integer                   :: MAXlpHalo! Max (over all PEs) lp halo size used
integer                   :: MAXmpHalo! Max (over all PEs) mp halo size used
integer                   :: comm
real*8                    :: time
character(64)             :: name     ! region name
character(80) :: string

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

! Print timing results to file named timing.mype
!ret = gptlpr (mype)

if(parallelBuild) then
!SMS$INSERT call GET_SMS_MPI_COMMUNICATOR(COMM)
!SMS$INSERT ret = gptlpr_summary(COMM)
else
  open(20,FILE='timing.summary')
  write(20,'(20x,"REGION NAME     MIN TIME     MAX TIME")')
  ret = gptlget_nregions (0, nregions)
  do n=0,nregions-1
    name = ' '
    ret = gptlget_regionname (0, n, name)
    ret = gptlget_wallclock (trim (name), 0, time)
    write(20,'(1x,A30,F13.3)') trim (name), time
  end do
  close(20)
endif

!ret = gptlprint_memusage ('Memory usage:')

print*,'IPE completed successfully'

end subroutine stop
