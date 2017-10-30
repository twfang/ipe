Program loadBalance
!This routine does LP load balancing for the IPE code.
!This code reads in timings (times) from a representative run of the IPE code parallelized on NLP or NMP points in the LP or MP direction only: In SMSnamelist "process_layout=NLP,1" or "process_layout=1,NMP".
!The representative run should be made with "barriersOn=f' in IPE.inp and "load_balance_on=f,f" in SMSnamelist.
! 10/24/17 IPEr504 theia_intel_parallel_170 CTIPINT time was almost indetical with and without barriers.
!Since the run is made parallelized with NLP(NMP) points in the LP(MP) direction, lpHaloSize(mpHaloSize)  must be <=1 which may require a smaller timestep.
!Currently lpHaloSize(mpHaloSize) is set in main/module_input_parameters.f90
!Currently flux_tube_solver times are used which are grepped from timing.* by:
!grep '   flux_tube_solver   ' timing.* | cut -c 1-7 --complement | tr -d ':' | awk '{$2=$4=""}1' | sort -h > times
!This code combines some of the LP(MP) points to create groups (bins) that have similar flux_tube_solver times.
!To create the bins this code uses the the max time over all processors ( maxval(times) ) as the target time (targetTime).
!Ideally all the bins would have a total time of targetTime for a minimum number of bins and perfect load balancing.
!In practice the times don't add up that way and it's generally good to allow some bins to have a total time (maxTimeAllowed) of slightly more than targetTime.
!The code has two input:
!       The number of points (NP) which is either NLP or NMP
!       "tolerance" which controls the size of maxtimeAllowed by : maxTimeAllowed = tolerance*targetTime
!tolerance must be >= 1.0. A good value for tolerance is 1.02. 
!If you want to decrease the number of bins, say from 61 to 60 increase tolerance say from 1.02 to 1.03.
!Conversely, if you want to increase the number of bins, say from 59 to 60 decrease tolerance.
!For parallel distribution the groups, not the lp points, will be apportioned among the processors.

implicit none
integer             :: NP                 !Number of lp(mp) points
integer             :: maxBins            !Maximum number of lp(mp) points in a group.
integer             :: i                  !Index for the  groups
integer             :: nt                 !Index for the time in a group
integer             :: d                  !Index for the direction
integer             :: n,ns,ne,nd         !Indexes for NP
integer             :: p                  !Processor number starting from 0
integer             :: Ncalls             !Number of calls for a processor
integer             :: istat              !Status for open
integer             :: Ngroups(2)         !Number of groups
integer             :: ng,ng1,ng2         !Number of groups indexes
integer             :: status             !Status of the read statement.
integer,allocatable :: layout (:,:)       !Starting lp(mp) in a group
integer,allocatable :: layouts(:)         !Starting lp(mp) in a group
real                :: minTime            !Minimum input time
real                :: meanTime           !Mean input time
real                :: maxTime            !Maximum input time
real                :: totalTime          !Sum of all the input times
real                :: targetTime         !The time we are shooting for
real                :: tolerance=1.02     !See explanation above.
real                :: maxdiff(2)         !Maximum difference between low and high for each direction
real                :: testTime           !The sum of times in a group for testing
real                :: spread1,spread2    !Spread of forward(backward) pass
real                :: maxTimeAllowed     !The maximum time allowed in a group.
real                :: dum1,dum2,dum3,dum4!Dummy values for reading file 'times'
real,allocatable    :: times(:)           !The run times which are read in from the 'ties' file
real,allocatable    :: sumTime (:,:)      !Running sum of the times in a group
real,allocatable    :: sumTimes(:)        !Temporary storage for sumTimes
character(LEN=13)   :: direction(2)       !The direction for the load balancing: forward or backward
character(LEN=72)   :: fname              !The name of the output file

direction(1) = 'forward pass '
direction(2) = 'backward pass'
print"('Enter number of points (NLP or NMP): ',$)"
read(*,*,iostat=status) NP
print"('Enter tolerance >= 1.0 or press Ctrl-D for the default value 1.02: ',$)"
read(*,*,iostat=status) tolerance
if(tolerance < 1.0) then
  print*,'tolerance must be >= 1.0, you entered ',tolerance
  stop
endif

maxBins = NP
allocate(times(NP))
allocate(layout(0:maxBins,2) , layouts(0:maxBins) , sumTime(maxBins,2) , sumTimes(maxBins))
open(88,file='times')
do n=1,NP
  read(88,*) p,Ncalls,times(n),dum1,dum2,dum3,dum4
  if(p /= n-1) then
    print*,'Error reading pe2s'
    print*,n,p,Ncalls,times(n),dum1,dum2,dum3,dum4
    stop
  endif
enddo

maxTime        = maxval(times)
minTime        = minval(times)
meanTime       = sum(times)/NP
totalTime      = sum(times)
targetTime     = maxtime
maxTimeAllowed = tolerance*targetTime
maxdiff        = 0.0
ns             = 1
ne             = NP
nd             = 1
do d=1,2 !First forward pass, then backward pass
  i=1
  sumTime(:,d) = 0.0
  layout (:,d) = 0
  N_LOOP:  do n=ns,ne,nd
    do nt=1,maxBins
      testTime = sumTime(i,d)+times(n)
      if( testTime>= targetTime) then !Above or equal to the target
        if (testTime-targetTime <= targetTime-sumTime(i,d).and.testTime <= maxTimeAllowed) then !Less above target than below so use this one
          if(i>maxBins) then
            print*,'i > maxBins',i,maxBins
            stop
          endif
          sumTime(i,d) = testTime
          maxDiff(d) = max(maxDiff(d),abs(sumTime(i,d)-targetTime))
          layout (i,d) = n
          i=i+1
          cycle N_LOOP
        else !more above target than below so use the last one
          layout (i,d) = n -1 + 2*(d/2)
          maxDiff(  d) = max(maxDiff(d),abs(sumTime(i,d)-targetTime))
          i=i+1
          cycle
        endif
      else ! Below target so continue
        if(n == ne) then !At the end of the times
          if(sumTime(i-1,d)+testTime-targetTime <= maxDiff(d)) then !Add testTime to the previous bin
            i = i-1
          endif
          if(d==1) then
            sumTime(i,d) = sumTime(i,d)+testTime
            maxDiff(  d) = max(maxDiff(d),abs(sumTime(i,d)-targetTime))
            layout (i,d) = n
          endif
          exit N_LOOP
        endif
      endif
      sumTime(i,d) = sumTime(i,d)+times(n)
      cycle N_LOOP
    enddo ! nt loop
    print*,'Exited nt loop without filling all the bins'
    print*,d,n,i
    stop
  enddo N_LOOP
  ns=NP
  ne=1
  nd=-1
  Ngroups(d) = i
! print*,' '
! print*,direction(d)
! print*,layout(1:Ngroups(d),d)
! print*,'Points =',Ngroups(d)
! print*,'maxDiff=', maxDiff(d)
! print*,maxval(abs(sumTime(1:Ngroups(d),d)-targetTime))
! print*,'Min    =', minval(sumTime(1:Ngroups(d),d))
! print*,'Max    =', maxval(sumTime(1:Ngroups(d),d))
! print*,'Spread =', maxval(sumTime(1:Ngroups(d),d)) - minval(sumTime(1:Ngroups(d),d))
! print*,sumTime(1:Ngroups(d),d)
enddo !d

if(maxDiff(1) < maxDIff(2)) then
  if(maxval(sumTime(1:Ngroups(1),1)) <= maxTimeAllowed) then
    d = 1
  else
    d = 2
  endif
elseif(maxDiff(1) == maxDIff(2)) then
  ng1=Ngroups(1)
  ng2=Ngroups(2)
  spread1 = maxval(sumTime(1:ng1,1)) - minval(sumTime(1:ng1,1))
  spread2 = maxval(sumTime(1:ng2,2)) - minval(sumTime(1:ng2,2))
  if(spread1 <= spread2) then
    d = 1
  else
    d = 2
  endif
else
  d = 2
endif

ng=Ngroups(d)
if(d==2) then
  layouts    = layout(:,2)
  layouts(0) = NP+1
  sumTimes   = sumTime(:,d)
  do i=1,ng
    layout (i,d) = layouts (ng-i)-1
    sumTime(i,d) = sumTimes(ng-i+1)
  enddo
else
  layout(0,d) = 0
endif

print*
print"(' NP              ',i10  )",NP
print"(' Min    time     ',f10.2)",minTime
print"(' Mean   time     ',f10.2)",meanTime
print"(' Max    time     ',f10.2)",maxTime
print"(' Total  time     ',f10.2)",totalTime
print"(' Target time     ',f10.2)",targetTime
print"(' tolerance       ',f10.5)",tolerance
print"(' Max time allowed',f10.2)",maxTimeAllowed
print*
print*,direction(d)
print"(' Number of Groups',i10  )", ng
print"(' Max Difference  ',f10.2)", maxDiff(d)
print"(' Min             ',f10.2)", minval(sumTime(1:ng,d))
print"(' Max             ',f10.2)", maxval(sumTime(1:ng,d))
print"(' Spread          ',f10.2)", maxval(sumTime(1:ng,d)) - minval(sumTime(1:ng,d))
print*
print*,'       bin     N range     time       spread'
do i=1,ng
  print"(i10,i9'-',i3,2f12.2)",i-1,layout(i-1,d)+1,layout(i,d),sumTime(i,d),sumTime(i,d)-targetTime
enddo

fname = './load_balance_groups'
OPEN (100,FILE=fname,STATUS='REPLACE',IOSTAT=istat)
write(100,*) layout(1:ng,d)

end Program loadBalance

