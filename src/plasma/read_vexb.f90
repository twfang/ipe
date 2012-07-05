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
!NOTE: 20120308: read the EXB drift values from old runs
      subroutine read_vexb ( utime,lp,mp )
        USE module_precision
        USE module_IPE_dimension,ONLY:NMP,NLP,IPDIM
        USE module_PLASMA,ONLY:VEXBup
        USE module_input_parameters,ONLY:start_time,stop_time
        USE module_IO, ONLY:filename,FORM_dum,STATUS_dum,luntmp1,luntmp2
        IMPLICIT NONE
        INTEGER (KIND=int_prec),  INTENT(IN)  :: utime    !universal time [sec]
        INTEGER (KIND=int_prec),INTENT(IN) :: lp,mp
!---local variables---
        INTEGER (KIND=int_prec) :: record_number_plasma_dum,utime_dum


luntmp1=6001
luntmp2=6002
if( utime==start_time ) then 
!ut_rec
          filename ='exb_ut_rec'
          FORM_dum ='formatted  ' 
          STATUS_dum ='old'
          CALL open_file ( filename, luntmp1, FORM_dum, STATUS_dum )

!plasma
          filename ='exb_plasma16'
          FORM_dum ='unformatted' 
          STATUS_dum ='old'
          CALL open_file ( filename, luntmp2, FORM_dum, STATUS_dum )
endif
      READ (UNIT=luntmp1,FMT=*) record_number_plasma_dum, utime_dum
print *,'check rec#',record_number_plasma_dum
print *,'check ut', utime_dum,' utime-172800=',(utime-172800)
if ( utime_dum /= (utime-172800) ) then
print *,'check ut', utime_dum,' utime-172800=',(utime-172800)
STOP
endif
!ExB
      READ (UNIT=luntmp2) VEXBup

!print *,'check EXB MAX',MAXVAL( vexbup(1:NMP,lp) ),MINVAL( vexbup(1:NMP,lp) )
print *,'check EXB lp=135',mp, vexbup(mp,135) ,utime

if( utime==stop_time ) then 
             CLOSE(UNIT=luntmp1)
             CLOSE(UNIT=luntmp2)
endif

END      subroutine read_vexb
