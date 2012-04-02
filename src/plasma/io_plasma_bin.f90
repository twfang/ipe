!CAUTION!!!!!!!the plasma i-o orders/file names have been changed !!!
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
      SUBROUTINE io_plasma_bin ( switch, utime )
      USE module_precision
      USE module_IO,ONLY: LUN_PLASMA1,LUN_PLASMA2,lun_min1,lun_min2,lun_ut,lun_ut2,record_number_plasma,lun_max1
      USE module_PLASMA,ONLY: plasma_3d,plasma_3d4n,VEXBup
      USE module_FIELD_LINE_GRID_MKS,ONLY: JMIN_IN,JMAX_IS
      USE module_IPE_dimension,ONLY: NMP0,NMP1,NLP,NPTS2D,ISPEC,ISPEV,NLP_all,IPDIM,ISPET
      USE module_input_parameters,ONLY:sw_debug,record_number_plasma_start
      USE module_physical_constants,ONLY:zero
      IMPLICIT NONE
!------------------------
      INTEGER (KIND=int_prec), INTENT(IN) :: switch !2:read; 1:write
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      REAL (KIND=real_prec),DIMENSION(:,:), ALLOCATABLE :: dumm  !(NPTS2D,NMP)
      INTEGER (KIND=int_prec) :: stat_alloc
      INTEGER (KIND=int_prec) :: jth,mp,lp,npts
      INTEGER (KIND=int_prec),pointer :: lun,in,is
      INTEGER (KIND=int_prec) :: n_read,n_read_min, utime_dum,record_number_plasma_dum
      INTEGER (KIND=int_prec),PARAMETER :: n_max=10000
!----------------------------------

IF ( switch<1.or.switch>2 ) THEN
  print *,'sub-io_plasma:!STOP! INVALID switch',switch
  STOP
END IF

if(sw_debug)  print *,'sub-io_plasma_bin: switch=',switch,' utime[sec]' ,utime 


!output time dependent plasma parameters
IF (.NOT.ALLOCATED(dumm) ) THEN
  ALLOCATE ( dumm(1:NPTS2D,NMP0:NMP1) &
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

!1:WRITE TO FILE
IF ( switch==1 ) THEN

record_number_plasma = record_number_plasma+1

j_loop1: DO jth=1,(ISPEC+3+ISPEV)

mp_loop1:do mp=NMP0,NMP1
 lp_loop1:do lp=1,nlp
!  i_loop:do i=
IN=>JMIN_IN(lp)
IS=>JMAX_IS(lp)
  npts = IS-IN+1 

IF ( jth<=ISPEC ) THEN
  dumm(IN:IS,mp) = plasma_3d(mp,lp)%N_m3(jth,        1:npts)
ELSE IF ( jth==(ISPEC+1) ) THEN
  dumm(IN:IS,mp) = plasma_3d(mp,lp)%Te_k(            1:npts)
ELSE IF ( jth<=(ISPEC+3) ) THEN
  dumm(IN:IS,mp) = plasma_3d(mp,lp)%Ti_k(jth-ISPEC-1,1:npts)
ELSE IF ( jth<=(ISPEC+3+ISPEV) ) THEN
  dumm(IN:IS,mp) = plasma_3d4n(IN:IS,mp)%V_ms1(jth-ISPEC-3)
ELSE 

END IF
! end do i_loop!i


!print *,'IN'
IF ( ASSOCIATED(IN) ) NULLIFY(in,is)
 end do lp_loop1!lp
end do mp_loop1!mp


LUN => LUN_PLASMA1(jth-1+lun_min1)
if(sw_debug) print *,'jth=',jth,' LUN=',LUN

      WRITE (UNIT=lun) dumm
if(sw_debug) print *,'!dbg! output dummy finished'

!print *,'lun',
IF ( ASSOCIATED(lun) ) NULLIFY(lun)
END DO j_loop1!jth


LUN => LUN_PLASMA1(lun_max1)
!ExB
      WRITE (UNIT=LUN) VEXBup
IF ( ASSOCIATED(LUN) ) NULLIFY(LUN)
!t if(sw_debug) print *,'!dbg! output VEXB finished'

!UT
      WRITE (UNIT=lun_ut,FMT=*) record_number_plasma, utime
if(sw_debug) & 
     &  print *,'LUN=',lun_ut,'!dbg! output UT  finished: utime=',utime,record_number_plasma




!2:READ FROM FILE
ELSE IF ( switch==2 ) THEN
print *,'sub-io_pl: start_uts=',utime

! array initialization
  DO mp=NMP0,NMP1
    DO lp=1,NLP_all
      DO npts=1,IPDIM
        plasma_3d(mp,lp)%N_m3(1:ISPEC,npts) = zero
        plasma_3d(mp,lp)%Te_k(        npts) = zero
        plasma_3d(mp,lp)%Ti_k(1:ISPET,npts) = zero
      END DO
    END DO
  END DO
!UT
      read_loop0: DO n_read=1,n_max !=10000
      READ (UNIT=lun_ut2,FMT=*) record_number_plasma_dum, utime_dum
print *,' record_# ', record_number_plasma_dum,' uts=', utime_dum
        IF (n_read==1) THEN
          n_read_min=record_number_plasma_dum
print *,'n_read=',n_read,'n_read_min=', n_read_min
        END IF
        IF (record_number_plasma_dum==record_number_plasma_start) THEN
           IF (utime_dum==utime) THEN
print *,'n_read=',n_read,' confirmed that ipe output exist at rec#=',record_number_plasma_dum,' at UT=', utime_dum
             CLOSE(UNIT=lun_ut2)
             EXIT  read_loop0
           ELSE !IF (utime_dum/=utime) THEN
print *,'n_read',n_read,'!STOP! INVALID start_time!',utime, record_number_plasma_dum, utime_dum
             STOP
           END IF

        ELSE IF (n_read==n_max.OR.record_number_plasma_dum>record_number_plasma_start) THEN
print *,'n_read',n_read,'!STOP! INVALID record number!',record_number_plasma_dum, utime_dum
          STOP
        END IF
      END DO read_loop0!: DO n_read=1,n_max          
        

j_loop2: DO jth=1,(ISPEC+3)  !t  +ISPEV)


LUN => LUN_PLASMA2(jth-1+lun_min2)
if(sw_debug) print *,'jth=',jth,' LUN2=',LUN

read_loop: DO n_read=n_read_min, record_number_plasma_start
if(jth==ISPEC+3) print *,'n_read=',n_read
!(1) UT
!      READ (UNIT=lun ) utime_dum
!if(sw_debug) 
!print *,'LUN=',lun,'!dbg! read UT  finished: jth=',jth,utime_dum
      READ (UNIT=lun ) dumm
!if(sw_debug) 
if (jth==ISPEC+3) print *,'!dbg! read dummy finished jth',jth


!      IF ( utime_dum==utime ) THEN
!        EXIT read_loop
!      ELSE
!        IF ( n_read==n_read_max ) THEN 
!          print *,'!STOP! INVALID reading plasma file',n_read,n_read_max,jth,LUN
!          STOP
!        END IF
!      END IF !      IF ( utime_dum==start_time ) THEN

    END DO read_loop !: DO n_read=1,n_read_max


mp_loop2:do mp=NMP0,NMP1
 lp_loop2:do lp=1,NLP
IN=>JMIN_IN(lp)
IS=>JMAX_IS(lp)
  npts = IS-IN+1 

IF ( jth<=ISPEC ) THEN
  plasma_3d(mp,lp)%N_m3(jth,         1:npts) = dumm(IN:IS,mp) 
ELSE IF ( jth==(ISPEC+1) ) THEN
  plasma_3d(mp,lp)%Te_k(             1:npts) = dumm(IN:IS,mp)
ELSE IF ( jth<=(ISPEC+3) ) THEN
  plasma_3d(mp,lp)%Ti_k(jth-ISPEC-1, 1:npts) = dumm(IN:IS,mp)
!t ELSE IF ( jth<=(ISPEC+3+ISPEV) ) THEN
!t   plasma_3d4n(IN:IS,mp)%V_ms1(jth-ISPEC-3  ) = dumm(IN:IS,mp)
END IF

!print *,'IN'
IF ( ASSOCIATED(IN) ) NULLIFY(in,is)
 end do lp_loop2!lp
end do mp_loop2!mp


print *,'closing lun',LUN
CLOSE(LUN)
IF ( ASSOCIATED(lun) ) NULLIFY(lun)
END DO j_loop2!jth

END IF !( switch==1 ) THEN

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

print *,'END sub-io_pl: sw=',switch,' uts=' ,utime 
      END SUBROUTINE io_plasma_bin
