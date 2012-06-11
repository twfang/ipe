!note:20120207: v36: used only activating the perp.transport gradually only during daytime...
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
!
      PROGRAM  test_plasma
      USE module_precision
      USE module_input_parameters,ONLY: read_input_parameters,start_time,stop_time,time_step,HPEQ_flip,ip_freq_msis,sw_output_plasma_grid ,sw_debug, sw_perp_transport, sw_rw_sw_perp_trans
      USE module_FIELD_LINE_GRID_MKS,ONLY: init_plasma_grid
      USE module_NEUTRAL_MKS,ONLY: neutral 
      USE module_PLASMA,ONLY: plasma
      USE module_ELDYN,ONLY: init_eldyn, eldyn
      USE module_IO,ONLY: open_output_files,output,close_files
      USE module_IPE_dimension,ONLY: NMP_all
      IMPLICIT NONE

      INTEGER (KIND=int_prec)   :: utime !universal time [sec]
      INTEGER(KIND=int_prec),parameter :: luntmp=300
      INTEGER(KIND=int_prec) :: istat,mp

      WRITE(*,*)" DATE: 08 September, 2011"
      WRITE(*,*)"********************************************"
      WRITE(*,*)"***      Copyright 2011 NAOMI MARUYAMA   ***"
      WRITE(*,*)"***      ALL RIGHTS RESERVED             ***"
      WRITE(*,*)"********************************************"
      WRITE(*,*)" LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model"
      WRITE(*,*)" DEVELOPER: Dr. Naomi Maruyama"
      WRITE(*,*)" CONTACT INFORMATION:"
      WRITE(*,*)" E-MAIL : Naomi.Maruyama@noaa.gov"
      WRITE(*,*)" PHONE  : 303-497-4857"
      WRITE(*,*)" ADDRESS: 325 Broadway, Boulder, CO 80305"
      WRITE(*,*)"                                            "

! set up input parameters
      CALL read_input_parameters ( )

! open Input/Output files
      CALL open_output_files ( )

! create allocatable arrays
      CALL allocate_arrays ( 0 )

! set up plasma grids by reading file
      CALL init_plasma_grid ( )


IF ( sw_output_plasma_grid ) THEN
  print *, 'sub-init_p: output plasma_grid'
  CALL output_plasma_grid ( )
END IF




! initialise the flux tubes from previous runs
      IF ( HPEQ_flip==0.0 ) THEN
        print *,'before CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time
        CALL io_plasma_bin ( 2, start_time )
        print *,'after CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time

      END IF
!20120215: CALL io_plasma_bin_readinit ( start_time )
!20120215: print *,'CALL io_plasma_bin_readinit finished!'


! initialization of electrodynamic module:
! read in E-field
      IF ( sw_perp_transport>=1 ) & 
     & CALL init_eldyn ( )

      time_loop: DO utime = start_time, stop_time, time_step

! updates auroral precipitation

! interplate from plasma to neutral grid: Nei,Tei,Vi,NHEAT, auroral heating?

!nm20110907:moved here because empirical Efield is needed for both neutral &plasma
      IF ( sw_perp_transport>=1 ) & 
     &   CALL eldyn ( utime )

! update neutral 3D structure: use MSIS/HWM to get the values in the flux tube grid
        IF ( MOD( (utime-start_time),ip_freq_msis)==0 ) THEN 
 IF ( sw_debug )  print *,'call MSIS',utime,start_time,ip_freq_msis,(utime-start_time),MOD( (utime-start_time),ip_freq_msis)
          CALL neutral ( utime )
        END IF

! interpolate from neutral to plasma grid:Tn,Un,[O,N2,O2],EHT(1,k), auroral heating?

! update plasma
        CALL plasma ( utime )


! update self-consistent electrodynamics
!t        CALL eldyn ( utime )

! output to a file
        CALL output ( utime )

      END DO  time_loop !: DO utime = start_time, stop_time, time_step

! DEallocate arrays
      CALL allocate_arrays ( 1 )

! close all open files
      CALL close_files ( )


!dbg20120509: no longer need 
!20120207 I have to output the sw_perp_transport to a file for the next run...
!IF ( sw_rw_sw_perp_trans ) THEN
!open(unit=luntmp, file='fort.300',status='unknown',form='formatted',iostat=istat)
!DO mp=1,NMP_all
!write(unit=luntmp, fmt='(2i3)') mp,sw_perp_transport(mp)
!print *,'mp=',mp,' sw_p',sw_perp_transport(mp)
!END DO
!close(unit=luntmp)
!END IF !( sw_tmp_sw_perp_trans ) THEN

END PROGRAM  test_plasma
