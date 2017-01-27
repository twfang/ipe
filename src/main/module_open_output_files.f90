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
      MODULE module_open_output_files
      USE module_precision
      USE module_IPE_dimension,ONLY: ISPEC,ISPEV
      IMPLICIT NONE

!nm20121003:module parameters are separated into module_io.f90!

      PRIVATE
      PUBLIC :: open_output_files

      CONTAINS
!---------------------------
        SUBROUTINE open_output_files ( )
        USE module_input_parameters,ONLY: NYEAR,NDAY,HPEQ_flip,sw_debug,sw_output_plasma_grid,record_number_plasma_start,sw_output_fort167
        USE module_IO,ONLY: &
&  filename,FORM_dum,STATUS_dum &
&, LUN_pgrid,PRUNIT,LUN_LOG &
&, LUN_FLIP1,LUN_FLIP2,LUN_FLIP3,LUN_FLIP4 &
&, LUN_PLASMA0, LUN_PLASMA1,LUN_PLASMA2, LUN_UT, LUN_UT2 &
&, lun_min1,lun_max1,lun_min2,lun_max2 &
&, record_number_plasma,luntmp1,luntmp2,luntmp3
        USE module_open_file,ONLY: open_file

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
END MODULE module_open_output_files

