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
      MODULE module_IO
      USE module_precision
      USE module_IPE_dimension,ONLY: ISPEC,ISPEV
      IMPLICIT NONE

! --- PRIVATE ---
!
! --- PUBLIC ---
      CHARACTER (LEN=100), PUBLIC :: filename
      CHARACTER (LEN=11),           PUBLIC :: FORM_dum
      CHARACTER (LEN=7),            PUBLIC :: STATUS_dum 
!--- unit numbers
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: LUN_pgrid=7
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: PRUNIT=8
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: LUN_LOG=9
      INTEGER (KIND=int_prec), PUBLIC :: LUN_flip1, LUN_flip2, LUN_flip3, LUN_flip4
      INTEGER (KIND=int_prec), PUBLIC :: LUN_PLASMA0, LUN_UT, LUN_UT2
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: lun_min1=100
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: lun_max1=lun_min1+ISPEC+3+ISPEV-1+1
      INTEGER (KIND=int_prec), DIMENSION(lun_min1:lun_max1),PUBLIC,TARGET :: LUN_PLASMA1 !WRITE
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: lun_min2=200
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: lun_max2=lun_min2+ISPEC+3-1
      INTEGER (KIND=int_prec), DIMENSION(lun_min2:lun_max2),PUBLIC,TARGET :: LUN_PLASMA2 !READ
      INTEGER (KIND=int_prec), PUBLIC :: record_number_plasma
!nm20120311
      INTEGER (KIND=int_prec), PUBLIC :: luntmp1,luntmp2,luntmp3
      PRIVATE
      PUBLIC :: open_output_files,output,close_files

      CONTAINS
!---------------------------
        SUBROUTINE open_output_files ( )
        USE module_input_parameters,ONLY: NYEAR,NDAY,HPEQ_flip,sw_debug,sw_output_plasma_grid,record_number_plasma_start,sw_output_fort167
        IMPLICIT NONE
        CHARACTER (LEN=100) :: string_tmp
        INTEGER (KIND=int_prec)::i
!
        IF ( NDAY < 100 ) THEN
          WRITE ( string_tmp, FMT="(i4,A2,i2)" ) NYEAR,'_0',NDAY
        ELSE         ! NDAY>=100
          WRITE ( string_tmp, FMT="(i4,A1,i3)" ) NYEAR,'_' ,NDAY
        END IF
        WRITE( UNIT=LUN_LOG, FMT=*) string_tmp

!--- unit=8
!nm20120303        filename ='FLIP_ERROR_FLAG_'//TRIM(string_tmp)//'.log'
        filename ='FLIP_ERR'
        FORM_dum ='formatted  ' 
        STATUS_dum ='unknown'
        CALL open_file ( filename, PRUNIT, FORM_dum, STATUS_dum )  

!--- unit=9
!nm20120303        filename ='logfile'//TRIM(string_tmp)//'.log'
        filename ='input_par'
        FORM_dum ='formatted  ' 
        STATUS_dum ='unknown'
        CALL open_file ( filename, LUN_LOG, FORM_dum, STATUS_dum )  

        IF ( sw_output_fort167 ) THEN
!--- unit=167
        LUN_FLIP1=167
        WRITE ( filename, FMT="('fort.',i3)" ) LUN_FLIP1  !'fort.167'
print *,'fort.167?', filename, LUN_FLIP1
        FORM_dum ='formatted  ' 
        STATUS_dum ='unknown'
        CALL open_file ( filename, LUN_FLIP1, FORM_dum, STATUS_dum ) 

!--- unit=168
        LUN_FLIP2=168
        WRITE ( filename, FMT="('fort.',i3)" ) LUN_FLIP2  !'fort.168'
print *,'fort.168?', filename, LUN_FLIP2
        CALL open_file ( filename, LUN_FLIP2, FORM_dum, STATUS_dum )

!--- unit=170
        LUN_FLIP3=170
        WRITE ( filename, FMT="('fort.',i3)" ) LUN_FLIP3  !'fort.170'
print *,'fort.170?', filename, LUN_FLIP3
        CALL open_file ( filename, LUN_FLIP3, FORM_dum, STATUS_dum )

!--- unit=171
        LUN_FLIP4=171
        WRITE ( filename, FMT="('fort.',i3)" ) LUN_FLIP4  !'fort.171'
print *,'fort.171?', filename, LUN_FLIP4
        CALL open_file ( filename, LUN_FLIP4, FORM_dum, STATUS_dum )
        END IF !( sw_debug ) THEN


IF ( sw_output_plasma_grid ) THEN
        LUN_PLASMA0=98
        FORM_dum = 'unformatted' 
        filename = 'plasma_grid'
        CALL open_file ( filename, LUN_PLASMA0, FORM_dum, STATUS_dum )
END IF !( sw_output_plasma_grid ) THEN

!nm20110923        LUN_PLASMA2=99
!nm20110923        filename = 'plasma_startup0'
!nm20110923        CALL open_file ( filename, LUN_PLASMA2, FORM_dum, STATUS_dum )

        LUN_UT=lun_min1-1 !=99
!nm20120303        filename ='ut_rec.log'
        filename ='ut_rec'
        FORM_dum ='formatted  ' 
        STATUS_dum ='unknown'
        CALL open_file ( filename, LUN_UT, FORM_dum, STATUS_dum ) 


!--- plasma output: unit=100~~115
        FORM_dum = 'unformatted' 
        DO i=lun_min1,lun_max1
           if(sw_debug) print *,'plasma: unit=',i        
           LUN_PLASMA1(i)=i
           IF ( (i-lun_min1) < 10 ) THEN
              WRITE( string_tmp, FMT="('0',i1)" )(i-lun_min1)
           ELSE IF ( (i-lun_min1) < 100 ) THEN
              WRITE( string_tmp, FMT="(i2)" )(i-lun_min1)
           END IF
           filename ='plasma'//TRIM(string_tmp)
           if(sw_debug) print *,(i-lun_min1),'filename',filename
           CALL open_file ( filename, LUN_PLASMA1(i), FORM_dum, STATUS_dum )
        END DO
        record_number_plasma = record_number_plasma_start - 1

        IF ( HPEQ_flip==0.0 ) THEN
!--- unit=1181 : input history file
!          LUN_PLASMA12=1181
!          STATUS_dum ='old'
!          filename = 'plasma_startup1'
!          CALL open_file ( filename, LUN_PLASMA12, FORM_dum, STATUS_dum )

          LUN_UT2=lun_min2-1 !=199
!nm20120303          filename ='startup_ut_rec.log'
          filename ='stup_ut_rec'
          FORM_dum ='formatted  ' 
          STATUS_dum ='old'
          CALL open_file ( filename, LUN_UT2, FORM_dum, STATUS_dum )

!--- unit=200~~215
          FORM_dum = 'unformatted' 
!note: reading Vi will be needed for neutral coupling
          DO i=lun_min2,lun_max2   !t +ISPEV)
             !if(sw_debug) 
             print *,'startup: unit=',i        
             LUN_PLASMA2(i)=i
             IF ( (i-lun_min2) < 10 ) THEN
                WRITE( string_tmp, FMT="('0',i1)" )(i-lun_min2)
             ELSE IF ( (i-lun_min2) < 100 ) THEN
                WRITE( string_tmp, FMT="(i2)" )(i-lun_min2)
             END IF
!nm20120303             filename ='startup'//TRIM(string_tmp)
             filename ='stup'//TRIM(string_tmp)
             !if(sw_debug) 
             print *,(i-lun_min2),'filename',filename
             CALL open_file ( filename, LUN_PLASMA2(i), FORM_dum, STATUS_dum )
          END DO
        END IF  ! ( HPEQ_flip==0.0 ) THEN

        END SUBROUTINE open_output_files
!---------------------------
        SUBROUTINE output ( utime )
        IMPLICIT NONE
!------------------------
        INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]

        WRITE(UNIT=PRUNIT,FMT="('uts=',i7)") utime

        END SUBROUTINE output
!---------------------------
        SUBROUTINE close_files ( )
        IMPLICIT NONE
        INTEGER (KIND=int_prec) :: i
!
        CLOSE(UNIT=PRUNIT)
        CLOSE(UNIT=LUN_UT)

DO i=lun_min1,lun_max1
        CLOSE( UNIT=LUN_PLASMA1(i) )
END DO

        END SUBROUTINE close_files
!---------------------------
END MODULE module_IO

!-------
        SUBROUTINE open_file ( filename_dum, UNIT_dum, FORM_dum, STATUS_dum )  
        USE module_precision
        USE module_IO,ONLY: LUN_LOG
        IMPLICIT NONE
        CHARACTER (LEN=*), INTENT(IN) :: filename_dum
        INTEGER (KIND=int_prec), INTENT(IN) :: UNIT_dum
        CHARACTER (LEN=*), INTENT(IN) :: FORM_dum
        CHARACTER (LEN=*), INTENT(IN) :: STATUS_dum
!---local
        LOGICAL :: flag
        INTEGER (KIND=int_prec) :: istat

IF (UNIT_dum>=1) THEN
  print *,'unit number=',UNIT_dum
ELSE
  print *,'!STOP! unit number not provided!!!'
  STOP
END IF

        INQUIRE ( UNIT=UNIT_dum, OPENED=flag ) 
WRITE( UNIT=LUN_LOG, FMT=*) flag
        IF ( .NOT. flag ) THEN
WRITE( UNIT=LUN_LOG, FMT=*) 'opening file:',filename_dum
          OPEN(UNIT=UNIT_dum,FILE=TRIM(filename_dum),STATUS=STATUS_dum,FORM=FORM_dum,IOSTAT=istat)
          IF ( istat /= 0 ) THEN
            WRITE( UNIT=LUN_LOG, FMT=*)'ERROR OPENING FILE',filename_dum
            STOP
          END IF
        END IF      

        END SUBROUTINE open_file
