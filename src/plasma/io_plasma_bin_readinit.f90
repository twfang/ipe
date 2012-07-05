!20120215: io_plasma_bin_readinit: special initialization to read the mp=1 file and distribute the values to all over mp... 
!date:Thu Sep 29 18:58:05 GMT 2011
! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!--------------------------------------------        
      SUBROUTINE io_plasma_bin_readinit ( utime )
      USE module_precision
      USE module_IO,ONLY: LUN_PLASMA1,LUN_PLASMA2,lun_min1,lun_min2,lun_ut,lun_ut2,record_number_plasma,lun_max1
      USE module_PLASMA,ONLY: plasma_3d,plasma_3d4n,VEXBup
      USE module_FIELD_LINE_GRID_MKS,ONLY: JMIN_IN,JMAX_IS, plasma_grid_3d,ISL,IBM,IGR,IQ,IGCOLAT,IGLON
      USE module_IPE_dimension,ONLY: NPTS2D,ISPEC,ISPEV,NLP,IPDIM,ISPET,NMP
      USE module_input_parameters,ONLY:sw_debug,record_number_plasma_start
      USE module_physical_constants,ONLY:zero,pi
      IMPLICIT NONE
!------------------------
!      INTEGER (KIND=int_prec), INTENT(IN) :: switch !2:read; 1:write
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      REAL (KIND=real_prec),DIMENSION(:,:), ALLOCATABLE :: dumm  !(NPTS2D,NMP)
      INTEGER (KIND=int_prec) :: stat_alloc
      INTEGER (KIND=int_prec) :: jth,mp,lp,npts
      INTEGER (KIND=int_prec),pointer :: lun,in,is
      INTEGER (KIND=int_prec) :: n_read,n_read_min, utime_dum,record_number_plasma_dum
      INTEGER (KIND=int_prec),PARAMETER :: n_max=10000
!20120215
      INTEGER (KIND=int_prec) :: lpj,midpoint,mpx,record_number_min
      INTEGER (KIND=int_prec),DIMENSION(NMP) :: record_number,flag
      INTEGER (KIND=int_prec),DIMENSION(ISPEC+3,NMP) :: flag2d
      REAL (KIND=real_prec),PARAMETER :: dlt=0.125 !within15min
      REAL (KIND=real_prec),PARAMETER :: utime_min=248306-86400+900  !sec: the begining of the final 24hrs
      REAL (KIND=real_prec),DIMENSION(NMP) :: ltime
      INTEGER (KIND=int_prec), parameter :: mp1=1
      REAL (KIND=real_prec) :: ltime_mp1
      INTEGER (KIND=int_prec), parameter :: n_read_min1=17 !utime=75506
      INTEGER (KIND=int_prec), parameter :: n_read_max=209!utime=248306
!----------------------------------

!IF ( switch/=3 ) THEN
!  print *,'sub-io_plasma_readinit:!STOP! INVALID switch',switch
!  STOP
!END IF
!
!print *,'sub-io_plasma_readinit: switch=',switch,' utime[sec]' ,utime 


!output time dependent plasma parameters
IF (.NOT.ALLOCATED(dumm) ) THEN
  ALLOCATE ( dumm(1:NPTS2D,NMP) &
     &,STAT=stat_alloc )         
      IF ( stat_alloc/=0 ) THEN
        print *,"sub-io_p:!STOP! ALLOCATION FAILD!:",stat_alloc
        STOP
      END IF
ELSE
STOP 'sub-io_p:!STOP! dumm has been allocated already???!!!'
END IF

! array initialization
dumm(:,:)=zero

!2:READ FROM FILE
!IF ( switch==3 ) THEN
print *,'start_time=',utime

! array initialization
  DO mp=1,NMP
    DO lp=1,NLP
      DO npts=1,IPDIM
        plasma_3d(mp,lp)%N_m3(1:ISPEC,npts) = zero
        plasma_3d(mp,lp)%Te_k(        npts) = zero
        plasma_3d(mp,lp)%Ti_k(1:ISPET,npts) = zero
      END DO
    END DO
  END DO


!jicamarca
  lpj=130
! calculate LT at midpoint
  midpoint = JMIN_IN(lpj) + ( JMAX_IS(lpj) - JMIN_IN(lpj) )/2
!calculate LT(mp) [hr] for utime @ midpoint
  do  mp=1,NMP
     ltime(mp)=REAL(utime)/3600. + plasma_grid_3d(midpoint,mp,IGLON)*180.0/pi/15.0
     IF ( ltime(mp) > 24.0 )  ltime(mp) = MOD(ltime(mp), 24.0)
print *,mp,'ltime',ltime(mp),' glon',(plasma_grid_3d(midpoint,mp,IGLON)*180.0/pi)
  end do

!UT
  record_number(1:NMP)=0
  flag(1:NMP)=0
  read_loop0: DO n_read=n_read_min1, n_read_max
     READ (UNIT=lun_ut2,FMT=*) record_number_plasma_dum, utime_dum



     print *,'n_read',n_read,' record_number_plasma_dum', record_number_plasma_dum,' utime_dum', utime_dum
     IF (utime_dum<utime_min) THEN
!      print *,'goto the next n_read! utime_dum',utime_dum,' utime_min=', utime_min
      CYCLE read_loop0
     END IF


!print *,'goto the next part! utime_dum=',utime_dum,' utime_min=', utime_min
!utime_dum
     ltime_mp1=REAL(utime_dum)/3600. + plasma_grid_3d(midpoint,mp1,IGLON)*180.0/pi/15.0
     IF ( ltime_mp1 > 24.0 )  ltime_mp1 = MOD(ltime_mp1, 24.0)
!print *, 'ltime_mp1', ltime_mp1


     mp_loop0:do mp=1,NMP
       IF ( (ltime_mp1-dlt) <= ltime(mp).AND. ltime(mp) <= (ltime_mp1+dlt) ) then
print *,'ltmp1-dlt',(ltime_mp1-dlt),' lt_mp',ltime(mp),' ltmp1+dlt',(ltime_mp1+dlt) 
         if ( record_number(mp)==0 ) then
           record_number(mp)=record_number_plasma_dum
           print *,' mp',mp,' glon=',plasma_grid_3d(midpoint,1,IGLON)*180.0/pi
           flag(mp)=1
           CYCLE read_loop0
!note: i hope that more than one record_number is matching,,,but error check just in case,,, 
         else if ( record_number(mp)>0 ) then
           print *,'!INVALID! record_number', mp,record_number(mp)
           STOP
         end if
       END IF
     end do mp_loop0


      END DO read_loop0!: DO n_read=1,n_max

!error check
      do mp=1,nmp
         if ( flag(mp)==0 ) then           
            print *,mp,flag(mp),'!STOP! matching record not found!!!'
            STOP
         end if
      end do


!determine record_number_min
record_number_min = MINVAL(record_number)
print *,'record_number_min',record_number_min 

flag2d(1:ISPEC+3 , 1:NMP)=0
j_loop2: DO jth=1,(ISPEC+3)  !t  +ISPEV)

LUN => LUN_PLASMA2(jth-1+lun_min2)
if(sw_debug) print *,'jth=',jth,' LUN2=',LUN

read_loop: DO n_read=n_read_min1, n_read_max
if(jth==ISPEC+3) print *,'n_read=',n_read


!(1) UT
!      READ (UNIT=lun ) utime_dum
!if(sw_debug) 
!print *,'LUN=',lun,'!dbg! read UT  finished: jth=',jth,utime_dum
      READ (UNIT=lun ) dumm
!if(sw_debug) 
if (jth==ISPEC+3) print *,'!dbg! read dummy finished jth',jth



! go to the next read if recordnumber < recordnumber_min
      if ( n_read<record_number_min ) then 
         CYCLE read_loop
      else

!check if recordnumber match?
         mp_loop1: do mp=1,nmp
            IF ( n_read==record_number(mp) ) THEN
               mpx=mp
               EXIT mp_loop1

            ELSE
               if ( mp==nmp ) then
                  print *,'matching record not found, goto the next read!'
                  CYCLE read_loop
               end if
            END IF
         end do mp_loop1

      end if !( n_read<record_number_min ) then 



!mp_loop2:do mp=1,NMP
      mp=mpx
      flag2d(jth,mp)=1
      print *,'assigning mp1 values at mp=',mp
      lp_loop2:do lp=1,NLP
         IN=>JMIN_IN(lp)
         IS=>JMAX_IS(lp)
         npts = IS-IN+1 

         IF ( jth<=ISPEC ) THEN
            plasma_3d(mp,lp)%N_m3(jth,         1:npts) = dumm(IN:IS,mp1) 
         ELSE IF ( jth==(ISPEC+1) ) THEN
            plasma_3d(mp,lp)%Te_k(             1:npts) = dumm(IN:IS,mp1)
         ELSE IF ( jth<=(ISPEC+3) ) THEN
            plasma_3d(mp,lp)%Ti_k(jth-ISPEC-1, 1:npts) = dumm(IN:IS,mp1)
!t ELSE IF ( jth<=(ISPEC+3+ISPEV) ) THEN
!t   plasma_3d4n(IN:IS,mp)%V_ms1(jth-ISPEC-3  ) = dumm(IN:IS,mp)
         END IF

!print *,'IN'
         IF ( ASSOCIATED(IN) ) NULLIFY(in,is)
      end do lp_loop2!lp
!end do mp_loop2!mp

   END DO read_loop !: DO n_read=1,n_read_max

print *,'closing lun',LUN
CLOSE(LUN)
IF ( ASSOCIATED(lun) ) NULLIFY(lun)
END DO j_loop2!jth

!END IF !( switch==1 ) THEN

IF ( ALLOCATED(dumm) ) THEN 
  DEALLOCATE ( dumm &
     &,STAT=stat_alloc )         
      IF ( stat_alloc/=0 ) THEN
        print *,"sub-io_p:!STOP! DEALLOCATION FAILD!:",stat_alloc
        STOP
      END IF
ELSE
STOP 'sub-io_p:!STOP! dumm has not been allocated???!!!'
END IF




!error check2
      do jth=1,ISPEC+3
      do mp=1,nmp
         if ( flag2d(jth,mp)==0 ) then           
            print *,jth,mp,flag2d(jth,mp),'!STOP!  values not assigned!!! ERROR!'
            STOP
         end if
      end do
      end do

print *,'END sub-io_plasma_bin_readinit '!,switch
      END SUBROUTINE io_plasma_bin_readinit
