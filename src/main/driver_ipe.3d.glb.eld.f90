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
!SMS$IGNORE BEGIN
      USE module_ELDYN,ONLY: init_eldyn, eldyn
!SMS$IGNORE END
      USE module_IO,ONLY: open_output_files,output,close_files
      USE module_IPE_dimension,ONLY: NMP_all,NLP_all
      IMPLICIT NONE
      include "gptl.inc"

      INTEGER (KIND=int_prec)   :: utime !universal time [sec]
      INTEGER(KIND=int_prec),parameter :: luntmp=300
      INTEGER(KIND=int_prec) :: istat,mp,ret

      call gptlprocess_namelist ('GPTLnamelist', 77, ret) 
      ret = gptlinitialize ()
      ret = gptlstart ('Total')

!SMS$CREATE_DECOMP(dh,<NMP_all,NLP_all>,<5,5>)

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
      ret = gptlstart ('read_input')
      CALL read_input_parameters ( )
      ret = gptlstop  ('read_input')

! open Input/Output files
      ret = gptlstart ('open_output_files')
      CALL open_output_files ( )
      ret = gptlstop  ('open_output_files')

! create allocatable arrays
      ret = gptlstart ('allocate_arrays0')
      CALL allocate_arrays ( 0 )
      ret = gptlstop ('allocate_arrays0')

! set up plasma grids by reading file
      ret = gptlstart ('init_plasma_grid')
      CALL init_plasma_grid ( )
      ret = gptlstop  ('init_plasma_grid')


IF ( sw_output_plasma_grid ) THEN
  ret = gptlstart ('output_plasma_grid')
  print *, 'sub-init_p: output plasma_grid'
  CALL output_plasma_grid ( )
  ret = gptlstop  ('output_plasma_grid')
END IF




! initialise the flux tubes from previous runs
      IF ( HPEQ_flip==0.0 ) THEN
        print *,'before CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time
        ret = gptlstart ('io_plasma_bin')
        CALL io_plasma_bin ( 2, start_time )
        ret = gptlstop  ('io_plasma_bin')
        print *,'after CALL io_plasma_bin finished! READ: start_time=', start_time,stop_time

      END IF
!20120215: CALL io_plasma_bin_readinit ( start_time )
!20120215: print *,'CALL io_plasma_bin_readinit finished!'


! initialization of electrodynamic module:
! read in E-field
      ret = gptlstart ('init_eldyn')
      IF ( sw_perp_transport>=1 ) &
     & CALL init_eldyn ( )
      ret = gptlstop  ('init_eldyn')

      ret = gptlstart ('time_loop')
      time_loop: DO utime = start_time, stop_time, time_step
      print*,'utime=',utime
! updates auroral precipitation

! interplate from plasma to neutral grid: Nei,Tei,Vi,NHEAT, auroral heating?

!nm20110907:moved here because empirical Efield is needed for both neutral &plasma
      ret = gptlstart ('eldyn')
      IF ( sw_perp_transport>=1 ) & 
     &   CALL eldyn ( utime )
      ret = gptlstop  ('eldyn')

! update neutral 3D structure: use MSIS/HWM to get the values in the flux tube grid
        IF ( MOD( (utime-start_time),ip_freq_msis)==0 ) THEN 
 IF ( sw_debug )  print *,'call MSIS',utime,start_time,ip_freq_msis,(utime-start_time),MOD( (utime-start_time),ip_freq_msis)
          ret = gptlstart ('neutral')
          CALL neutral ( utime )
          ret = gptlstop  ('neutral')
        END IF

! interpolate from neutral to plasma grid:Tn,Un,[O,N2,O2],EHT(1,k), auroral heating?

! update plasma
        ret = gptlstart ('plasma')
        CALL plasma ( utime )
        ret = gptlstop  ('plasma')


! update self-consistent electrodynamics
!t        CALL eldyn ( utime )

! output to a file
        ret = gptlstart ('output')
        CALL output ( utime )
        ret = gptlstop  ('output')

      END DO  time_loop !: DO utime = start_time, stop_time, time_step
      ret = gptlstop  ('time_loop')

! DEallocate arrays
      ret = gptlstart ('allocate_arrays1')
      CALL allocate_arrays ( 1 )
      ret = gptlstop  ('allocate_arrays1')

! close all open files
      ret = gptlstart ('close_files')
      CALL close_files ( )
      ret = gptlstop  ('close_files')


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

      ret = gptlstop  ('Total')
      call stop

END PROGRAM  test_plasma
