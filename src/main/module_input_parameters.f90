!note:20120207: v36: used only activating the perp.transport gradually...
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
      MODULE module_input_parameters
      USE module_precision
      USE module_IPE_dimension,ONLY: NLP,NMP,NPTS2D
      IMPLICIT NONE

!--- IPE wide run parameters
      INTEGER (KIND=int_prec), PUBLIC   :: utime           !UT[sec] IPE internal time management
      INTEGER (KIND=int_prec), PUBLIC   :: nTimeStep=1     !internal number of time steps
      INTEGER (KIND=int_prec), PUBLIC   :: start_time      !=0  !UT[sec]
      INTEGER (KIND=int_prec), PUBLIC   :: stop_time       !=60 !UT[sec]
      INTEGER (KIND=int_prec), PUBLIC   :: time_step       !=60 ![sec]
      INTEGER (KIND=int_prec), PUBLIC   :: nprocs=1        !Number of processors
      INTEGER (KIND=int_prec), PUBLIC   :: mype=0          !Processor number
      INTEGER (KIND=int_prec), PUBLIC   :: lps,lpe,mps,mpe !Per processor start and stop indexes for lp,mp
      INTEGER (KIND=int_prec), PUBLIC   :: lpHaloSize=99   !lp halo size (big number=NOP for serial)
      INTEGER (KIND=int_prec), PUBLIC   :: mpHaloSize=99   !mp halo size (big number=NOP for serial)
      INTEGER (KIND=int_prec), PUBLIC   :: MaxLpHaloUsed=0 !Max lp halo size used for the entire run
      INTEGER (KIND=int_prec), PUBLIC   :: MaxMpHaloUsed=0 !Max mp halo size used for the entire run

      REAL (KIND=real_prec), PUBLIC :: F107D   !.. Daily F10.7
      REAL (KIND=real_prec), PUBLIC :: F107AV  !.. 81 day average F10.7
!
      INTEGER (KIND=int_prec), PUBLIC :: NYEAR ! year
      INTEGER (KIND=int_prec), PUBLIC :: NDAY  ! day number

      INTEGER (KIND=int_prec), PUBLIC :: internalTimeLoopMax=1  ![times]internal time loop: default 1
      INTEGER (KIND=int_prec), PUBLIC :: ip_freq_output=900  ![sec] must be multiple of time_step: default 15m
      INTEGER (KIND=int_prec), PUBLIC :: ip_freq_msis=180    !frequency[sec] to call MSIS/HWM: default 3m
      INTEGER (KIND=int_prec), PUBLIC :: ip_freq_plasma=60   !frequency[sec] to call plasma: default 1m
      INTEGER (KIND=int_prec), PUBLIC :: ip_freq_eldyn=180   !frequency[sec] to call eldyn: default 3m(for quiet climatology),60s for storm
      LOGICAL                , PUBLIC :: parallelBuild=.false.

!--- FLIP specific input parameters
      REAL (KIND=real_prec), PUBLIC :: DTMIN_flip  !.. Minimum time step allowed (&=10 secs?)
      INTEGER (KIND=int_prec),PUBLIC :: sw_INNO  !.. switch to turn on FLIP NO calculation if <0
      REAL (KIND=real_prec), PUBLIC :: FPAS_flip   !.. Pitch angle scattering fraction
      REAL (KIND=real_prec), PUBLIC :: HPEQ_flip   !.. Sets initial equatorial [H+][cm-3] if positive
      REAL (KIND=real_prec), PUBLIC :: HEPRAT_flip !.. Initial He+/H+ ratio (.01 to 1.0)
      REAL (KIND=real_prec), PUBLIC :: COLFAC_flip !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)

      INTEGER (KIND=int_prec),PUBLIC :: sw_TEI    !.. switches Te/Ti solutions ON if > 0
      INTEGER (KIND=int_prec),PUBLIC :: sw_OHPLS  !.. switches O+/H+ solutions ON if > 0

      INTEGER (KIND=int_prec),PUBLIC :: sw_IHEPLS !.. switches He+ diffusive solutions on if > 0
      INTEGER (KIND=int_prec),PUBLIC :: sw_INPLS  !.. switches N+ diffusive solutions on if > 0
      INTEGER (KIND=int_prec),PUBLIC :: sw_wind_flip  !.. switch for neutral wind input to FLIP: 1:ON; 0:ZERO wind
      INTEGER (KIND=int_prec),PUBLIC :: sw_depleted_flip  !.. switch for depleted flux tube in FLIP: 1:ON; 0:OFF 
      INTEGER (KIND=int_prec),PUBLIC :: start_time_depleted !.. time UT to start to deplete the flux tube
      INTEGER,PUBLIC :: sw_neutral_heating_flip  !.. switch for neutral heating calculation in FLIP: 1:ON; 0:OFF
      REAL (KIND=real_prec), PUBLIC :: init_Te_max !.. max Te[K] in the initial profile
      INTEGER, PUBLIC :: sw_DEBUG_flip           !.. switch to turn on debug writes:0=off; 1=on for solver
      INTEGER, PUBLIC :: sw_ERSTOP_flip          !.. switch to turn on STOP in SUB-WRITE_EFLAG
!dbg20120301: N+ BAND solver issue
      LOGICAL, PUBLIC :: sw_LCE  !local chemical equilibruium below ht_LCE[km]
      REAL (KIND=real_prec), PUBLIC :: ht_LCE !.. max ht[km] for LCE
      REAL (KIND=real_prec), PUBLIC :: ZLBNP_inp !.. ZLBNP
!dbg20120304:
      REAL (KIND=real_prec), PUBLIC :: FNFAC_flip !.. FNFAC in RSPRIM.FOR
!dbg20121129
      LOGICAL, PUBLIC :: sw_optw_flip !=F  !chemical routine is called before He+ solution for too inflated o+ density due to exb drift
!dbg20121130
      LOGICAL, PUBLIC :: sw_init_guess_flip !=F  !this might help in finding a solution for convergence error???
      INTEGER (KIND=int_prec), PUBLIC :: dt_init_guess_flip=60 !max DT for changing init_guess 
!dbg20121130
      REAL (KIND=real_prec), PUBLIC :: ZLBDY_flip=120.  !Lower boundary altitude

!--- MSIS/HWM specific input parameters
      REAL (KIND=real_prec), DIMENSION(7), PUBLIC :: AP   ! magnetic index(daily)
!.. MSIS: or when sw(9)=-1. :                 
!           - array containing:                                         
!             (1) daily ap                                              
!             (2) 3 hr ap index for current time                        
!             (3) 3 hr ap index for 3 hrs before current time           
!             (4) 3 hr ap index for 6 hrs before current time           
!             (5) 3 hr ap index for 9 hrs before current time           
!             (6) average of eight 3 hr ap indicies from 12 to 33 hrs pr
!                    to current time                                    
!             (7) average of eight 3 hr ap indicies from 36 to 57 hrs pr
!                    to current time 
!.. HWM:        ap - two element array with                                    
!             ap(1) = magnetic index(daily) (use 4 in lower atmos.)     
!             ap(2)=current 3hr ap index (used only when sw(9)=-1.) 
!
!--- ELDYN specific input parameters
      REAL (KIND=real_prec), PUBLIC :: kp_eld   ! geomagnetic index
!--- all the SWITCHes either integer or logical or character
      LOGICAL, PUBLIC :: sw_debug
      LOGICAL, PUBLIC :: sw_debug_mpi
      LOGICAL, PUBLIC :: sw_output_fort167 =.false.
      LOGICAL, PUBLIC :: sw_output_wind    =.false. !unit=6000,6001
      LOGICAL, PUBLIC :: barriersOn=.false. !true means turn on barriers.
      INTEGER(KIND=int_prec), PUBLIC :: peFort167=0 !default mype=0
      INTEGER(KIND=int_prec), PUBLIC :: mpfort167 = 10
      INTEGER(KIND=int_prec), PUBLIC :: lpfort167 = 14
      INTEGER(KIND=int_prec), DIMENSION(2), PUBLIC :: iout
      INTEGER(KIND=int_prec), PUBLIC :: mpstop=80
      INTEGER(KIND=int_prec), PUBLIC :: sw_neutral=3    
!0: WAM debug: use which ever ESMF fields are coming across for debugging purpose
!1: WAM default science mode: specify ESMF fields you wish to use
!2: GT
!3: MSIS(default)
!4: read in files
      LOGICAL, dimension(7), PUBLIC :: swNeuPar!     =.false. !f:OFF (from MSIS); t:ON (from WAM)
!determines which neutral parameters to derive from WAM only when sw_neutral=0/1? 
!1:tn; 2:un1(east); 3:un2(north); 4:un3(up); 5:[O]; 6:[O2]; 7:[N2]
      LOGICAL, PUBLIC :: swEsmfTime =.false.
      INTEGER(KIND=int_prec), PUBLIC :: sw_eldyn
!0:self-consistent eldyn solver; 1:WACCM efield ;2:  ;3: read in external efield
      INTEGER(KIND=int_prec), PUBLIC :: sw_3DJ !1:calculate 3D currents, je_3d
      INTEGER(KIND=int_prec), PUBLIC :: sw_pcp        !0:heelis; 1:weimer
      INTEGER(KIND=int_prec), PUBLIC :: sw_grid       !0:APEX; 1:FLIP
! if sw_grid=1 
!dbg20120304: nolonger used
!nm20120304      REAL (KIND=real_prec), PUBLIC :: PCO_flip  
!nm20120304      REAL (KIND=real_prec), PUBLIC :: BLON_flip 
      LOGICAL, PUBLIC :: sw_output_plasma_grid
!JFM  LOGICAL, PUBLIC :: sw_rw_sw_perp_trans
      LOGICAL, PUBLIC :: sw_dbg_perp_trans
      INTEGER(KIND=int_prec), PUBLIC :: sw_perp_transport 
!0:WITHOUT perpendicular transport
!1:THETA only transport included
!2:both THETA&PHI:transport included, NH/SH flux tubes are moving together with the same ExB drift
!3:both THETA&PHI:transport included, NH/SH flux tubes are moving separately with different ExB drift
! if sw_perp_tr=>1
      INTEGER (KIND=int_prec), PUBLIC :: lpmin_perp_trans !=15 :mlatN=78deg???
      INTEGER (KIND=int_prec), PUBLIC :: lpmax_perp_trans !=151:mlatN=5.64deg
      INTEGER (KIND=int_prec), PUBLIC :: sw_th_or_r
      INTEGER (KIND=int_prec), PUBLIC :: record_number_plasma_start
      INTEGER (KIND=int_prec), PUBLIC :: sw_record_number
!nm20160329: used only when HPEQ_flip=0.5
      INTEGER (KIND=int_prec), PUBLIC :: ut_start_perp_trans=0
      INTEGER (KIND=int_prec), PUBLIC :: duration !used when sw_record_n=1
      INTEGER (KIND=int_prec), PUBLIC :: sw_exb_up
! (0) self consistent electrodynamics
! (1) WACCM E empirical model
! (2) GIP empirical model
! (3) SUPIM empirical model
      INTEGER(KIND=int_prec), PUBLIC :: sw_para_transport 
!0:WITHOUT parallel transport (no calling to flux tube solver)
!1:parallel transport included
      INTEGER(KIND=int_prec), PUBLIC :: sw_ksi
!0: ksi_factor=1.0---no compressional effect/no adiabatic heating
!1: ksi_factor from richards thesis---including compressional effect/adiabatic heating ,,,used until 20120314 with interpolate_ft.v16.f90
!2:20120330: new way of calculating the ksi_factor with interpolate_ft.v17.f90
      INTEGER(KIND=int_prec), PUBLIC :: sw_divv
!0: div * V//=0
!1: div * V// included in the Te/i solver
!dbg20120313 
      REAL(KIND=real_prec), PUBLIC :: fac_BM
      LOGICAL,PUBLIC :: check_halo_on,compare_var_on,exact_parallel_sum,load_balance_on(2),set_process_layout !Read from SMSnamelist
      INTEGER(KIND=int_prec), PUBLIC :: compare_var_ntasks_1,compare_var_ntasks_2,load_balance_method(2),load_balance_size(2),process_layout(2) !Read from SMSnamelist
!
! MPI communicator to be passed to SMS
      integer, PUBLIC :: my_comm
!---
      NAMELIST/IPEDIMS/NLP,NMP,NPTS2D 
      NAMELIST/NMIPE/start_time &
     &,stop_time &
     &,time_step &
     &,F107D   &
     &,F107AV  &
     &,NYEAR  &
     &,NDAY   &
     &,internalTimeLoopMax &
     &,ip_freq_eldyn &
     &,ip_freq_output &
     &,ip_freq_msis &
     &,ip_freq_plasma
      NAMELIST/NMFLIP/DTMIN_flip  & 
     &,sw_INNO   & 
     &,FPAS_flip   & 
     &,HPEQ_flip   & 
     &,HEPRAT_flip & 
     &,COLFAC_flip & 
     &,sw_TEI &
     &,sw_OHPLS &
     &,sw_IHEPLS &
     &,sw_INPLS  &
     &,sw_wind_flip &
     &,sw_depleted_flip &
     &,start_time_depleted &
     &,sw_neutral_heating_flip &
     &,init_Te_max &
     &,sw_DEBUG_flip &
     &,sw_ERSTOP_flip &
     &,sw_LCE &
     &,ht_LCE &
     &,ZLBNP_inp &
     &,FNFAC_flip &
     &,sw_optw_flip &
     &,sw_init_guess_flip &
     &,dt_init_guess_flip &
     &,ZLBDY_flip 
      NAMELIST/NMMSIS/AP  &
     &,kp_eld
      NAMELIST/NMSWITCH/&
           &  sw_neutral     &
           &, swNeuPar       &
           &, swEsmfTime     &
           &, sw_eldyn     &
           &, sw_3DJ         &
           &, sw_pcp         &
           &, sw_grid        &
           &, sw_output_plasma_grid        &
!JFM       &, sw_rw_sw_perp_trans &
           &, sw_dbg_perp_trans &
           &, sw_perp_transport &
           &, lpmin_perp_trans &
           &, lpmax_perp_trans &
           &, sw_th_or_r &
           &, sw_exb_up &
           &, sw_para_transport &
           &, sw_ksi &
           &, sw_divv &
           &, mpstop  &
           &, sw_debug       &
           &, sw_debug_mpi   &
           &, sw_output_fort167   &
           &, sw_output_wind   &
           &, mpfort167   &
           &, lpfort167   &
           &, peFort167   &
           &, record_number_plasma_start   &
           &, sw_record_number   &
           &, ut_start_perp_trans   &
           &, duration   &
           &, fac_BM   &
           &, iout     &
           &, barriersOn
      NAMELIST/smsnamelist/ check_halo_on,compare_var_ntasks_1,compare_var_ntasks_2,compare_var_on,exact_parallel_sum,load_balance_method,load_balance_on,load_balance_size,process_layout,set_process_layout

!nm20120304           &, PCO_flip       &
!nm20120304           &, BLON_flip      &


      PRIVATE
      PUBLIC :: read_input_parameters


      CONTAINS
!---------------------------
! initialise plasma grids
        SUBROUTINE read_input_parameters ( )
        USE module_IPE_dimension,ONLY: NLP,NMP,NPTS2D
        IMPLICIT NONE
!---------
!MPI requirement 
!SMS$INSERT      include "mpif.h"
!---
        INTEGER(KIND=int_prec),PARAMETER :: LUN_nmlt=1
        INTEGER(KIND=int_prec),PARAMETER :: LUN_SMS_nmlt=88
        CHARACTER(LEN=*),PARAMETER :: INPTNMLT='IPE.inp'
        INTEGER(KIND=int_prec) :: IOST_OP=0
        INTEGER(KIND=int_prec) :: IOST_RD=0
        INTEGER (KIND=int_prec), PARAMETER :: LUN_LOG0=10  !output4input parameters only
        CHARACTER (LEN=*), PARAMETER :: filename='logfile_input_params.log'
        INTEGER (KIND=int_prec) :: istat        

!SMS$IGNORE BEGIN
        OPEN(LUN_nmlt,FILE=INPTNMLT,ERR=222,IOSTAT=IOST_OP,STATUS='OLD')
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=IPEDIMS  ,ERR=222,IOSTAT=IOST_RD)
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=NMIPE    ,ERR=222,IOSTAT=IOST_RD)
!SMS$IGNORE END

!SMS$INSERT lpHaloSize=2
!SMS$INSERT mpHaloSize=1
!
!set up MPI communicator for SMS
!(1) when NEMS is not used, pass MPI_COMM_WORLD into SET_COMMUNICATOR()
!SMS$INSERT        my_comm=MPI_COMM_WORLD
!(2) when NEMS is used, my_comm=mpiCommunicator has been assigned already in sub-myIPE_Init
!        print *, 'sub-read_input_para:my_comm=', my_comm
!SMS$SET_COMMUNICATOR( my_comm )
!
!nm20160608 sms debug
!SMS$CREATE_DECOMP(dh,<NLP,NMP>,<lpHaloSize,mpHaloSize>: <NONPERIODIC, PERIODIC>)
!!!SMS$CREATE_DECOMP(dh,<NLP,NMP>,<lpHaloSize,mpHaloSize>: <PERIODIC, PERIODIC>)

!SMS$SERIAL BEGIN
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=NMFLIP   ,ERR=222,IOSTAT=IOST_RD)
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=NMMSIS   ,ERR=222,IOSTAT=IOST_RD)
        REWIND LUN_nmlt
        READ(LUN_nmlt,NML=NMSWITCH ,ERR=222,IOSTAT=IOST_RD)

        OPEN(UNIT=LUN_LOG0,FILE=filename,STATUS='unknown',FORM='formatted',IOSTAT=istat)
        IF ( istat /= 0 ) THEN
          WRITE( UNIT=6, FMT=*)'ERROR OPENING FILE',filename
          STOP
        END IF
        WRITE(UNIT=LUN_LOG0, NML=NMIPE)
        WRITE(UNIT=LUN_LOG0, NML=NMFLIP)
        WRITE(UNIT=LUN_LOG0, NML=NMMSIS)
        WRITE(UNIT=LUN_LOG0, NML=NMSWITCH)

        WRITE(UNIT=LUN_LOG0,FMT=*)'NMP=',NMP,' NLP=',NLP,' NPTS2D=',NPTS2D

        WRITE(UNIT=LUN_LOG0,FMT=*)'real_prec=',real_prec,' int_prec=',int_prec

        CLOSE(LUN_LOG0)
        process_layout = 1
!SMS$INSERT        OPEN  (LUN_SMS_nmlt,FILE='SMSnamelist',ERR=222,IOSTAT=IOST_OP,STATUS='OLD')
!SMS$INSERT        REWIND LUN_SMS_nmlt
!SMS$INSERT        READ  (LUN_SMS_nmlt,NML=smsnamelist  ,ERR=222,IOSTAT=IOST_RD)
!SMS$INSERT        CLOSE (LUN_SMS_nmlt)
print*,'process_layout',process_layout
!SMS$SERIAL END
!SMS$INSERT        if(.not.set_process_layout) then
!SMS$INSERT          print*,'In SMSnamelist set_process_layout must be true.',set_process_layout
!SMS$INSERT          stop
!SMS$INSERT        endif
        CLOSE(LUN_nmlt)
222     IF ( IOST_OP /= 0 ) THEN
          WRITE(UNIT=LUN_nmlt, FMT=*) "OPEN NAMELIST FAILED!", IOST_OP
          STOP
        ELSEIF ( IOST_RD /= 0 ) THEN
          WRITE(UNIT=LUN_nmlt, FMT=*) "READ NAMELIST FAILED!", IOST_RD
          STOP
        ENDIF

stop_time=start_time+duration

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

!SMS$insert call NNT_NPROCS(nprocs)
!SMS$insert call NNT_ME    (mype  )
!SMS$TO_LOCAL(dh:<1,lps:lbound>,<1,NLP:ubound>) BEGIN
lps = 1
lpe = NLP
!SMS$TO_LOCAL END
!SMS$TO_LOCAL(dh:<2,mps:lbound>,<2,NMP:ubound>) BEGIN
mps = 1
mpe = NMP
!SMS$TO_LOCAL END
print *,'finished reading namelist:',filename
print *,' '
print"(' NLP:                 ',I6)",NLP
print"(' NMP:                 ',I6)",NMP
print"(' mpstop:              ',I6)",mpstop
print"(' stop_time            ',I6)",stop_time
print"(' Number of Processors:',I6)",nprocs
print"(' lpHaloSize:          ',I6)",lpHaloSize
print"(' mpHaloSize:          ',I6)",mpHaloSize
print *,' '
print *,' '

!dbg20120509        IF ( sw_rw_sw_perp_trans )  CALL setup_sw_perp_transport ()
!note:20120207: v36: used only activating the perp.transport gradually...


!dbg20160711
!SMS$IGNORE begin
print*,mype,'sub-read_input: swNeuPar',swNeuPar
!SMS$IGNORE end
 
        END SUBROUTINE read_input_parameters

END MODULE module_input_parameters
