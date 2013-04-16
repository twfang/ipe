!---------------------------------
!purpose: test the GT for IPE
!---------------------------------

PROGRAM  test_GT

USE module_precision

USE moduleAmplitude  ! has amplitude type
USE moduleTidalPhase ! has tidalPhase type
USE moduleSwitches  ! has the switches type
USE moduleDriverDebug   ! for calling debugging routines


USE module_input_parameters, ONLY: read_input_parameters, &
                             start_time, stop_time, time_step, HPEQ_flip, sw_neutral  
                             ! add nday  **********************

!------------------------------------------------------------------------
! Don't need this now, b/c I am reading this grid
! from other files, will need it later once we couple to IPE completely
!------------------------------------------------------------------------
!USE module_FIELD_LINE_GRID_MKS, ONLY: init_plasma_grid  
! lrm20120531
!USE module_FIELD_LINE_GRID_MKS, ONLY: JMIN_IN, JMAX_IS  ! For testing purposes only *****


USE module_NEUTRAL_MKS, ONLY: neutral 

USE moduleTHERMOSPHERE, ONLY : GT_thermosphere_init, low_lat_efield, &
                               Foster, GT_thermosphere, &
                               calculate_magnetic_parameters_using_apex

USE modSizeFluxTube, ONLY : NPTS, NMP, NLP  ! sizes of flux tube grid


USE moduleInterfaceThermo2Iono, ONLY : INTERFACE__thermosphere_to_FIXED_GEO, &
                            INTERFACE__FIXED_GEO_to_IONOSPHERE

USE modSizeFixedGridIono, ONLY : nFixedGridIonoHeights, &
                                 nFixedGridIonoLats, &
                                 nFixedGridIonoLons

USE modSizeFixedGridThermo, ONLY : nFixedGridThermoHeights

! sizes of thermosphere pressure grid
USE modSizeThermo, ONLY : GT_ht_dim, GT_lat_dim, GT_lon_dim


USE moduleInterfaceIono2Thermo, ONLY : INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO, &
                                       readIPEtoGeoGrid, INTERFACE__FIXED_GRID_to_THERMO

      
!t    USE module_PLASMA,ONLY: plasma
!t    USE module_ELDYN,ONLY: init_eldyn, eldyn
!t    USE module_IO,ONLY: open_output_files,output,close_files

IMPLICIT NONE

INTEGER (KIND=int_prec)   :: utime ! universal time [sec]
      
      
! Variables for GT_thermosphere -------------------  NOW THIS IN A MODULE
! Parameters
!integer, parameter :: GT_ht_dim = 15
!integer, parameter :: GT_lat_dim = 91
!integer, parameter :: GT_lon_dim = 20


! Parameters for IPE fixed grid  - in modSizeFixedGridIono module
! nFixedGridIonoHeights 


! check status of fortran open
INTEGER :: FileOpenStatus


!------------------------------------------------------
! File unit number for checking startup IPE values :
!------------------------------------------------------
! THIS IS NOT BEING USED RIGHT NOW,  NEED TO WRITE TO THIS FILE
! INSTEAD OF PRINTING TO THE SCREEN ****************************** lrm20120830
INTEGER, parameter :: unitCheckStartIPE = 19
! Debug the Thermospheric values???
LOGICAL, parameter :: debugStartupIPE = .FALSE.
! THIS IS NOT BEING USED RIGHT NOW,  NEED TO WRITE TO THIS FILE
! INSTEAD OF PRINTING TO THE SCREEN ****************************** lrm20120830
CHARACTER(LEN=*), PARAMETER :: debugStartIPEFileName = 'CheckStartIPE.dat'


! ONLY FOR TESTING................
INTEGER(KIND=int_prec), POINTER :: IN, IS  !  take this out later lrm20120531


!------------------------------------------------------
! File unit number for checking thermosphere values :
!------------------------------------------------------
INTEGER, parameter :: unitCheckThermo =7700

!------------------------------------
! Debug the Thermospheric values???
!------------------------------------
LOGICAL, parameter :: debugThermo = .FALSE.
CHARACTER(LEN=*), PARAMETER :: debugThermoFileName = 'CheckGTGIP.dat'


!--------------------------------------------------------------------
! File unit number for checking Thermospheric interpolation values :
!--------------------------------------------------------------------
INTEGER, parameter :: unitCheckThermoInterp = 7800
!---------------------------------------------------
! Write out the Thermospheric interpolated values??
!---------------------------------------------------
LOGICAL, parameter :: debugThermoInterp = .FALSE.
CHARACTER(LEN=*), PARAMETER :: debugThermoInterpFileName = 'interpOut.dat'
INTEGER, parameter :: unitCheckThermoInterpBefore = 7600  ! for writing out thermosphere values
INTEGER, parameter :: unitCheckThermoInterpAfter = 7500  ! for writing out thermosphere values
INTEGER, parameter :: unitCheckThermoInterpFluxT = 7400  ! for writing out thermosphere values
INTEGER, parameter :: unitCheckPressureHeightThermo = 8800

!-----------------------------------------------------------------
! File unit number for checking Ionospheric interpolation values :
!-----------------------------------------------------------------
! THIS IS NOT BEING USED RIGHT NOW,  NEED TO WRITE TO THIS FILE
INTEGER, parameter :: unitCheckIonoInterpBefore = 9800
INTEGER, parameter :: unitCheckIonoInterpAfter = 9900

INTEGER, parameter :: numIonoVars = 10
!INTEGER :: unitCheckIonoInterp(numIonoVars)
CHARACTER(LEN=30) :: ionoVarName(numIonoVars)

!-----------------------------------------------------------------------------
! File unit number for checking Ionospheric values to the pressure grid
! and file unit number for getting the average height for each pressure level
!-----------------------------------------------------------------------------
INTEGER, parameter :: unitCheckPressureInterp = 8900
INTEGER, parameter :: unitCheckPressureHeight = 7820
INTEGER, parameter :: unitPressureGridHeight = 7810

!------------------------------------------------
! Write out the Ionospheric interpolated values??
!------------------------------------------------
LOGICAL, parameter :: debugIonoInterp = .FALSE.
! THIS IS NOT BEING USED RIGHT NOW,  NEED TO WRITE TO THIS FILE
! INSTEAD OF PRINTING TO THE SCREEN ****************************** lrm20120830
!CHARACTER(LEN=*), PARAMETER :: debugIonoInterpFileName = 'interpIonoOut.dat'


!---------------------------------------------------------------
! Write out the Ionospheric pressure grid interpolated values??
!---------------------------------------------------------------
LOGICAL, parameter :: debugIonoFixedtoPressure = .FALSE.

!----------------------------------------
! File unit number for IPE startup files
!----------------------------------------
INTEGER, parameter :: unitStartUp = 211

!----------------------------------
! Output variables from GT model
!----------------------------------
REAL(kind=8) :: wind_southwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: wind_eastwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: wvz_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)


REAL(kind=8) :: rmt_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: Temperature_K_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: ht_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: O_density_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: O2_density_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: N2_density_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)



INTEGER, parameter :: houghSize = 181
      
REAL*8 :: hough11(houghSize)
REAL*8 :: hough22(houghSize)
REAL*8 :: hough23(houghSize)
REAL*8 :: hough24(houghSize)
REAL*8 :: hough25(houghSize)


REAL*8 :: solar_declination_angle_radians
REAL*8 :: universal_time_seconds  
REAL*8 :: Universal_Time_hours
REAL*8 :: Start_time_UT_hours   !,Stop_time_UT_hours

!-------------------------
! For GT_thermosphere  :
!-------------------------
INTEGER :: foster_level
REAL*8 :: foster_power
REAL*8 :: emaps(21,GT_lon_dim,7) , cmaps(21,GT_lon_dim,7)
REAL*8 :: magnetic_latitude_degrees(91,20)
REAL*8 :: magnetic_longitude_degrees(91,20)

REAL*8 :: B_magnitude_apex_nT(GT_lat_dim,GT_lon_dim)
REAL*8 :: B_dip_angle_apex_degrees(GT_lat_dim,GT_lon_dim)      
REAL*8 :: profile(GT_ht_dim,21)

REAL*8 :: qion3d(GT_ht_dim,GT_lat_dim,GT_lon_dim)


REAL*8 :: exns(2,45,GT_lon_dim) , eyns(2,45,GT_lon_dim) , ezns(2,45,GT_lon_dim)

INTEGER :: nn=0 , nnloop=0 , &
           nn_composition_counter=0 , nn_smoothing_counter=0

INTEGER :: ii, jj, kk, ll, mm

!INTEGER :: iitime = 1 ! FOR TESTING gt->fixed grid interpolation

INTEGER :: bigLoop = 1 ! for counting how many times we've been in outer loop
INTEGER :: littleLoop = 1 ! for counting how many times we've been in inner loop

INTEGER (KIND=int_prec) :: gtLoopTime

INTEGER :: idump_gt  ! dump out GT?


!----------------------------------------------
! Ionospheric Variables in thermospheric grid
!----------------------------------------------
REAL*8 :: Oplus_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL*8 :: Hplus_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)  ! NOT USED YET lrm20120508
REAL*8 :: N2plus_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim) ! NOT USED YET lrm20120508
REAL*8 :: Te_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)             ! NOT USED YET lrm20120508
REAL*8 :: Nplus_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)  ! NOT USED YET lrm20120508
REAL*8 :: Ne_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL*8 :: NOplus_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL*8 :: O2plus_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL*8 :: Ti_Oplus_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)


!-------------------------------------------------------------------------------
! Ionospheric Variables from IPE model Ne_, NOplus_, O2plus_, Ti_Oplus
!-------------------------------------------------------------------------------
REAL(kind=8) :: Oplus_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: Hplus_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: Heplus_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: Nplus_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: NOplus_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: O2plus_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: N2plus_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: Oplus2D_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: Oplus2P_density_from_IPE(NPTS, NMP)
REAL(kind=8) :: Te_from_IPE(NPTS, NMP)  ! electron temperature
REAL(kind=8) :: Ti_from_IPE(NPTS, NMP)  ! ion temperature


!--------------------------------------------------
! sum of all the densities : 
! O+, H+, He+, N+, NO+, O2+, N2+, O+(2D), O+(2P)
!--------------------------------------------------
REAL(kind=8) :: Ne_density_from_IPE(NPTS, NMP)  


!-------------------------------------------------------------
! Ionospheric Variables interpolated to fixed grid (mid_lat)
!-------------------------------------------------------------
REAL(kind=8) :: Oplus_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: Hplus_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: Nplus_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: NOplus_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: O2plus_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: N2plus_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: Ne_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: Te_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)
REAL(kind=8) :: Ti_high_res_fixed(nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons)


!-------------------------------------
! Switches for Ionospheric parameters
!-------------------------------------
logical ::  GIP_switches(GT_lon_dim)


REAL*8 :: Foster_potential(23, GT_lon_dim, 7)
REAL*8 :: Foster_Efield_amplification
REAL*8 :: plvu(49) , zonal(49)

INTEGER :: number_of_GT_time_steps_in_24_hours


!--------------------------------------------------
! For interpolating from Thermosphere to IPE grid
!--------------------------------------------------
CHARACTER*10, parameter :: thermospheric_model_name = 'CMAT2'
real*8 :: therm_model_geo_long_deg(GT_lon_dim)
real*8 :: therm_model_geo_lat_deg(GT_lat_dim)
real*8 :: Altitude_m_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)


real*8 :: Vn_Southwards_ms1_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: Vn_Eastwards_ms1_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: Vn_Upwards_ms1_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: Tn_K_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: O_density_m3_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: O2_density_m3_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: N2_density_m3_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: elx_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: ely_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)
real*8 :: qion3d_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)


!-------------------------------- 
! Fixed grid neutral atmosphere
!--------------------------------                                              
REAL(kind=8) :: O_density_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: O2_density_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: N2_density_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: V_East_FixedHeight(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: V_South_FixedHeight(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: V_Upward_FixedHeight(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: tts_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: qion3d_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: elx_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: ely_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)


INTEGER :: fileNPTS, fileNMP, fileNLP  ! read from IPE grid file, but not used
REAL (KIND=8), DIMENSION(1:NMP+1) :: mlon_rad
REAL (KIND=8), DIMENSION(NPTS) :: ipeGL  ! read from IPE grid file, but not used

INTEGER :: inGIP1d(NLP), isGIP1d(NLP)   ! 1d in, read from the IPE grid file
! in - n means northpole, is - s means southpole
INTEGER :: inGIP(NMP, NLP), isGIP(NMP, NLP)  ! 2d inGIP, isGIP grids to match up w/ old GIP subroutine call

REAL(kind=8) :: pz_plasma_1d(NPTS)  ! 1d in, read from the IPE grid file
REAL(kind=8) :: glat_plasma_3d(npts,nmp), glond_plasma_3d(npts,nmp), pz_plasma_3d(npts,nmp)

!---------------------
! For IPE grid read :
!---------------------
!INTEGER :: iwrite_plasma_interface not used lrm20121115


REAL(kind=8) :: TN_plasma_input_3d(npts, nmp), O_plasma_input_3d(npts, nmp), &
                O2_plasma_input_3d(npts, nmp), N2_plasma_input_3d(npts, nmp)

REAL(kind=8) :: V_east_plasma(npts, nmp), &
                V_south_plasma(npts, nmp), &
                V_upward_plasma(npts, nmp)

!------------------------------------------------------
! Output from fixed_geo 
! lat, lon, ht points to be used for 
! interpolation from thermo fixed grid to ipe grid
!------------------------------------------------------
INTEGER :: ilon1_3d_fixed_ht(npts,nmp), ilon2_3d_fixed_ht(npts,nmp)
INTEGER :: ilat1_3d_fixed_ht(npts,nmp), ilat2_3d_fixed_ht(npts,nmp)
INTEGER :: ispecial_3d_fixed_ht(npts,nmp)
INTEGER :: ihl_3d_fixed_ht(npts,nmp), ihu_3d_fixed_ht(npts,nmp)

! must always be true to start with
LOGICAL :: isFirstCallFixedHeight = .TRUE.

!---------------------
! Physical Parameters
!---------------------
REAL(kind=8), PARAMETER  :: PI = 3.141592654
REAL(kind=8), PARAMETER  :: DTR = PI/180.0 ! Degrees to Radians conversion
REAL(kind=8), PARAMETER  :: R0 = 6.370E06  ! in meters

!-------------------------------
! For converting apex grid data
!-------------------------------
REAL(kind=8) :: gr(NPTS,nmp)
REAL(kind=8) :: gcol(NPTS,nmp)
REAL(kind=8) :: glon(NPTS,nmp)

!----------------------------------------------------------------------
! Debug results from INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE ?
!----------------------------------------------------------------------
INTEGER, parameter :: unitFixedGeo = 15
CHARACTER(LEN=*), PARAMETER :: debugFixedGeoFileName = 'CheckFixedGeo.dat'
LOGICAL :: debugFixedGeo = .FALSE.



!------------------------------------
! File with the IPE grid coordinates
!------------------------------------
CHARACTER(LEN=*), PARAMETER :: plasmaGridFileName = 'plasma_grid.new.ascii'
INTEGER, PARAMETER :: unitGridIPE = 16

!------------------------------------------------------
! File for checking IPE grid values :
!------------------------------------------------------
INTEGER, parameter :: unitCheckIPE = 17
CHARACTER(LEN=*), PARAMETER :: debugIPEFileName = 'checkIPEgrid.txt'
!-------------------------------
! Debug the IPE grid values???
!-------------------------------
LOGICAL, parameter :: debugGridIPE = .FALSE.




!----------------------
! Files & Directories
!----------------------
CHARACTER(len=140) :: geoToMagFileName
CHARACTER(len=200) :: GT_input_dataset , GT_output_dataset
character(len=140) :: ionoStartDir
character(len=140) :: staticFileDir
character(len=140) :: debugDir
CHARACTER(len=140) :: giptogeoFileName


! should be inputs???
INTEGER :: nnstop , nnstrt
INTEGER :: nday
INTEGER :: GT_timestep_in_seconds
REAL*8  :: f107

!-------------
! inputs
!-------------
TYPE(amplitudeType)  :: amplitude   ! replaces ampl11, ...
TYPE(tidalPhaseType) :: tidalPhase  ! replaces lt11, ...

INTEGER :: smoothingFrequency
INTEGER :: neutralCompositionFrequency

TYPE(switchesType) :: switches   ! replaces sw_, .......

!---------------------------
! IPE Startup files type
!---------------------------
TYPE startUpType
  INTEGER :: unit
  CHARACTER(30) :: filename
  CHARACTER(7) :: speciesName
  REAL(kind=8) :: species(NPTS, NMP)
END TYPE startUpType

integer, parameter :: numIonoStart = 11


!! TO+ = TH+ = Ti  - just read above Ti & use for both TO+, TH+


type (startUpType) :: startUpFiles(numIonoStart)

!-------------------------
! For timing the code :
!-------------------------
REAL :: startTime, endTime
REAL :: timeFluxGridtoFixedGrid = 0.
REAL :: timeFixedGridtoPressureGrid = 0.
REAL :: timePressureGridtoFixedGrid = 0.
REAL :: timeFixedGridtoFluxGrid = 0.


!-----------------------------------
! Namelist for input parameters :
!-----------------------------------
NAMELIST /gtipeINPUTS/GT_input_dataset, GT_output_dataset, &
                      ionoStartDir, staticFileDir, debugDir, &
                      nday, f107, GT_timestep_in_seconds, &
                      amplitude, tidalPhase, switches, &
                      smoothingFrequency, neutralCompositionFrequency


! BEGIN CODE ====================================================================================

debugDir = TRIM(debugDir)


! Set GIP switches

!------------------------
! Set up input parameters
!------------------------
CALL read_input_parameters ()

!-------------------------
! Open Input/Output files
!-------------------------
!t   CALL open_output_files ()

!---------------------------------------
! Set up plasma grids by reading file
!---------------------------------------
! CALL init_plasma_grid ( )  ! Will need this later 


! initialise the flux tubes from previous runs
!t      IF ( HPEQ_flip==0.0 ) &
!t     &   CALL init_flux_tubes ( )

! initialization of electrodynamic module:
!t     uCALL init_eldyn ( )



!-------------------------------
! OPEN THE GT-IPE NAMELIST FILE
!-------------------------------
OPEN (UNIT=16, FILE='nameListGTIPE.txt', STATUS='OLD')
read (16, NML=gtipeINPUTS)

print *,'driver_ipe_gt.3d : amplitude = ',amplitude
print *,'driver_ipe_gt.3d : tidalPhase = ',tidalPhase

print *,'driver_ipe_gt : GT_input_dataset = ', GT_input_dataset

startUpFiles%speciesName = (/ "Oplus", "Hplus", "Heplus", "Nplus", "NOplus", "O2plus", &
                              "N2plus", "0plus2D", "0plus2P", "startTe", "startTi" /)


startUpFiles%filename = (/ "plasma00.cmb.ascii", "plasma01.cmb.ascii", "plasma02.cmb.ascii", &
                           "plasma03.cmb.ascii", "plasma04.cmb.ascii", "plasma05.cmb.ascii", &
		           "plasma06.cmb.ascii", "plasma07.cmb.ascii", "plasma08.cmb.ascii", &
			   "plasma09.cmb.ascii", "plasma10.cmb.ascii"/)


			   

!! TO+ = TH+ = Ti  - just read above Ti & use for both TO+, TH+
		   			
			   
!--------------------------------------
! Open all start up ionospheric files
!--------------------------------------
!------------------------------------------------
! Open IPE species files
! One for each species, with all time steps in it
!------------------------------------------------
do ii = 1, numIonoStart

   startUpFiles(ii)%unit = unitStartUp + (ii - 1)
   print *,'startUpFiles(ii)%unit = ', startUpFiles(ii)%unit
   print *,TRIM(GT_input_dataset)//TRIM(startUpFiles(ii)%filename)
   
   OPEN(UNIT=startUpFiles(ii)%unit, FILE=TRIM(ionoStartDir)//TRIM(startUpFiles(ii)%filename), &
              STATUS='old',FORM='formatted')

enddo ! ii




! ***************************************************************

IF ( sw_neutral == 'GT' ) then

    nnstrt = 1  ! ********************
    nnstop = 15  ! ******************** NEED TO CHANGE THIS - LRM
      

    ! These are in the namelist :
    !GT_input_dataset
    !!GT_output_dataset 
    !GT_output_dataset 
    !staticFileDir
    
    !--------------------------------------
    ! Define name of the geo to mag file
    !--------------------------------------
    giptogeoFileName = TRIM(staticFileDir)//'/GIP_Fixed_GEO_grid_lowres'

    !-----------------------------------------
    ! How often to smooth (in # of time steps)
    !-----------------------------------------
    !smoothingFrequency = 5

    !neutralCompositionFrequency = 1

    !---------------------------------------------------------------------------
    ! These will need to be in a namelist (or some sort of input)
    ! Right now they are hard wired to values in the variable declarations
    !---------------------------------------------------------------------------
    !F    --     sw_electro - electrodynamics
    !F    --     sw_External_model_provides_NO_N4S_densities
    !F    --     sw_External_model_provides_low_lat_E_fields
    !T    --     sw_input_Auroral_production_is_single_overall_rate


    !------------------------
    ! Initialize arrays :
    !------------------------
    solar_declination_angle_radians = 0.
    hough11 = 0.
    hough22 = 0.
    hough23 = 0.
    hough24 = 0.
    hough25 = 0.              ! OUTPUT

    !------------------------
    ! Initialize arrays :
    !------------------------
    wind_southwards_ms1_FROM_GT = 0.  ! read from startup file             ! OUTPUT
    wind_eastwards_ms1_FROM_GT = 0.   ! read from startup file             ! OUTPUT
    wvz_FROM_GT = 0.  ! upward nuetral wind                                ! OUTPUT
    rmt_FROM_GT = 0. ! composition                                         ! OUTPUT
    Temperature_K_FROM_GT = 0. ! nuetral temperature in 3D pressure coords ! OUTPUT
    ht_FROM_GT = 0. ! height at each grid point in 3D                      ! OUTPUT

    ! # of GT time steps in one day, based on 60 second time step ************
    number_of_GT_time_steps_in_24_hours =  1440 ! hardwired for now *****************

    !------------------------------
    ! read in tiros/foster data
    !------------------------------
    OPEN(21,FILE=TRIM(staticFileDir)//'ionprof',STATUS='old')
    READ(21,99001) emaps
    READ(21,99001) cmaps
    CLOSE(21)


    OPEN(23,FILE=TRIM(staticFileDir)//'prof2',STATUS='old')
    READ(23,99001) profile
    CLOSE(23)


    !------------------------------------------------
    ! Read the foster potential from the holt file
    !------------------------------------------------
    OPEN(25,FILE=TRIM(staticFileDir)//'holt',STATUS='old')
    READ(25,99001) Foster_potential
    CLOSE(25)
 
99001 FORMAT (1x,6E13.6)

    !-------------------------------------------------------
    ! from tucan_time.f90 :
    ! Open the hough input files to read the constants
    !-------------------------------------------------------
    OPEN(33,FILE=TRIM(staticFileDir)//'hough',STATUS='old')
    OPEN(36,FILE=TRIM(staticFileDir)//'hough11',STATUS='old')
    OPEN(40,FILE=TRIM(staticFileDir)//'hough25',STATUS='old')


    !------------------
    ! Get the IPE grid 
    ! *eventually get rid of this, IPE will provide the grid in the future
    !------------------
    OPEN (unitGridIPE, FILE=TRIM(staticFileDir)//TRIM(plasmaGridFileName), STATUS='OLD')
    READ (UNIT=unitGridIPE,FMT="(i10)") fileNMP
    READ (UNIT=unitGridIPE,FMT="(i10)") fileNLP
    READ (UNIT=unitGridIPE,FMT="(i10)") fileNPTS
    READ (UNIT=unitGridIPE,FMT="(20i10)") inGIP1d(1:NLP)
    READ (UNIT=unitGridIPE,FMT="(20i10)") isGIP1d(1:NLP)
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") mlon_rad( 1: NMP+1 ) !rad
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") Pz_plasma_1d(1:NPTS) ! meters
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") ipeGL(1:NPTS) ! radians
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") GLAT_plasma_3d(1:NPTS, 1:NMP)  ! radians
    READ (UNIT=unitGridIPE,FMT="(20E12.4)") GLONd_plasma_3d(1:NPTS, 1:NMP) ! radians
    CLOSE (unitGridIPE)

    print *,'driver_ipe_gt.3d : Done reading the IPE grid.'

    if (debugGridIPE) then 

        CALL checkGridIPE(debugDir, debugIPEFileName, unitCheckIPE, NMP, NLP, &
                          inGIP1d, isGIP1d,  mlon_rad)

        STOP

    endif ! debugGridIPE
    
   !--------------------------------------------------------------
   ! Open files for writing out interpreted variables b/f & after
   !--------------------------------------------------------------
   If (debugIonoInterp) then 
       OPEN (unitCheckIonoInterpBefore, FILE=TRIM(debugDir)//TRIM('OplusBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter, FILE=TRIM(debugDir)//TRIM('OplusAfter.txt'), STATUS='REPLACE')

       OPEN (unitCheckIonoInterpBefore+1, FILE=TRIM(debugDir)//TRIM('HplusBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter+1, FILE=TRIM(debugDir)//TRIM('HplusAfter.txt'), STATUS='REPLACE')

       OPEN (unitCheckIonoInterpBefore+2, FILE=TRIM(debugDir)//TRIM('NplusBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter+2, FILE=TRIM(debugDir)//TRIM('NplusAfter.txt'), STATUS='REPLACE')

       OPEN (unitCheckIonoInterpBefore+3, FILE=TRIM(debugDir)//TRIM('NOplusBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter+3, FILE=TRIM(debugDir)//TRIM('NOplusAfter.txt'), STATUS='REPLACE')

       OPEN (unitCheckIonoInterpBefore+4, FILE=TRIM(debugDir)//TRIM('O2plusBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter+4, FILE=TRIM(debugDir)//TRIM('O2plusAfter.txt'), STATUS='REPLACE')

       OPEN (unitCheckIonoInterpBefore+5, FILE=TRIM(debugDir)//TRIM('N2plusBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter+5, FILE=TRIM(debugDir)//TRIM('N2plusAfter.txt'), STATUS='REPLACE')

       OPEN (unitCheckIonoInterpBefore+6, FILE=TRIM(debugDir)//TRIM('NeBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter+6, FILE=TRIM(debugDir)//TRIM('NeAfter.txt'), STATUS='REPLACE')

       OPEN (unitCheckIonoInterpBefore+7, FILE=TRIM(debugDir)//TRIM('TeBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter+7, FILE=TRIM(debugDir)//TRIM('TeAfter.txt'), STATUS='REPLACE')

       OPEN (unitCheckIonoInterpBefore+8, FILE=TRIM(debugDir)//TRIM('TiBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckIonoInterpAfter+8, FILE=TRIM(debugDir)//TRIM('TiAfter.txt'), STATUS='REPLACE')

   end if !  (debugIonoInterp)

   !--------------------------------------------------------------
   ! Open files for writing out interpreted variables b/f & after
   !--------------------------------------------------------------
   If (debugIonoFixedtoPressure) then 
       OPEN (unitCheckPressureInterp, FILE=TRIM(debugDir)//TRIM('OplusPressureGrid.txt'), STATUS='REPLACE')
       OPEN (unitCheckPressureInterp+1, FILE=TRIM(debugDir)//TRIM('HplusPressureGrid.txt'), STATUS='REPLACE')
       OPEN (unitCheckPressureInterp+2, FILE=TRIM(debugDir)//TRIM('NplusPressureGrid.txt'), STATUS='REPLACE')
       OPEN (unitCheckPressureInterp+3, FILE=TRIM(debugDir)//TRIM('NOplusPressureGrid.txt'), STATUS='REPLACE')
       OPEN (unitCheckPressureInterp+4, FILE=TRIM(debugDir)//TRIM('O2plusPressureGrid.txt'), STATUS='REPLACE')
       OPEN (unitCheckPressureInterp+5, FILE=TRIM(debugDir)//TRIM('N2plusPressureGrid.txt'), STATUS='REPLACE')
       OPEN (unitCheckPressureInterp+6, FILE=TRIM(debugDir)//TRIM('NePressureGrid.txt'), STATUS='REPLACE')
       OPEN (unitCheckPressureInterp+7, FILE=TRIM(debugDir)//TRIM('TePressureGrid.txt'), STATUS='REPLACE')
       OPEN (unitCheckPressureInterp+8, FILE=TRIM(debugDir)//TRIM('TiPressureGrid.txt'), STATUS='REPLACE')

       OPEN (unitCheckPressureHeight, FILE=TRIM(debugDir)//TRIM('averagePressureHeightGrid.txt'), STATUS='REPLACE')


  end if ! (debugIonoFixedtoPressure) 

   !--------------------------------------------------------------
   ! Open files for writing out interpreted variables b/f & after
   !--------------------------------------------------------------
   If (debugThermoInterp) then 
       OPEN (unitCheckThermoInterpBefore, FILE=TRIM(debugDir)//TRIM('SouthWindBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter, FILE=TRIM(debugDir)//TRIM('SouthWindAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT, FILE=TRIM(debugDir)//TRIM('SouthWindFluxTube.txt'), STATUS='REPLACE')

       OPEN (unitCheckThermoInterpBefore+1, FILE=TRIM(debugDir)//TRIM('EastWindBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+1, FILE=TRIM(debugDir)//TRIM('EastWindAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+1, FILE=TRIM(debugDir)//TRIM('EastWindFluxTube.txt'), STATUS='REPLACE')

       OPEN (unitCheckThermoInterpBefore+2, FILE=TRIM(debugDir)//TRIM('WVZBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+2, FILE=TRIM(debugDir)//TRIM('WVZAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+2, FILE=TRIM(debugDir)//TRIM('WVZFluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+3, FILE=TRIM(debugDir)//TRIM('RMTBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+3, FILE=TRIM(debugDir)//TRIM('RMTAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+3, FILE=TRIM(debugDir)//TRIM('RMTFluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+4, FILE=TRIM(debugDir)//TRIM('TemperatureBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+4, FILE=TRIM(debugDir)//TRIM('TemperatureAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+4, FILE=TRIM(debugDir)//TRIM('TemperatureFluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+5, FILE=TRIM(debugDir)//TRIM('OBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+5, FILE=TRIM(debugDir)//TRIM('OAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+5, FILE=TRIM(debugDir)//TRIM('OFluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+6, FILE=TRIM(debugDir)//TRIM('O2Before.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+6, FILE=TRIM(debugDir)//TRIM('O2After.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+6, FILE=TRIM(debugDir)//TRIM('O2FluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+7, FILE=TRIM(debugDir)//TRIM('N2Before.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+7, FILE=TRIM(debugDir)//TRIM('N2After.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+7, FILE=TRIM(debugDir)//TRIM('N2FluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+8, FILE=TRIM(debugDir)//TRIM('QionBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+8, FILE=TRIM(debugDir)//TRIM('QionAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+8, FILE=TRIM(debugDir)//TRIM('QionFluxTube.txt'), STATUS='REPLACE')

       OPEN (unitCheckPressureHeightThermo, FILE=TRIM(debugDir)//TRIM('averagePressureHeightGridfromThermo.txt'), STATUS='REPLACE')

       ! Write all heights out to a file 
       OPEN (unitPressureGridHeight, FILE=TRIM(debugDir)//TRIM('PressureHeightGrid.txt'), STATUS='REPLACE')

   end if !  (debugThermoInterp)

    !------------------------------------------------------------------------------
    ! Need to change 1d inGIP, isGIP, to 2d to match up with the old GIP subroutine 
    !------------------------------------------------------------------------------
    DO ii = 1, NLP
       inGIP(:,ii) = inGIP1d(ii)
       isGIP(:,ii) = isGIP1d(ii)
    ENDDO ! ii


    !------------------------------------------------------------------------------
    ! Convert heights in meters to height in kilometers
    !------------------------------------------------------------------------------
    DO ii = 1, NPTS
       Pz_plasma_3d(ii,:) = Pz_plasma_1d(ii)/1000.
    ENDDO ! ii


    !--------------------------------------------------------
    ! Convert co-latitude in radians to latitude in degrees
    !--------------------------------------------------------
    GLAT_plasma_3d = 90. - (GLAT_plasma_3d/DTR)


    !--------------------------------------------------------
    ! Convert longitude in radians to degrees
    !--------------------------------------------------------
    GLONd_plasma_3d = GLONd_plasma_3d/DTR


    !----------------------
    ! Debug the IPE grid
    !----------------------
    if (debugGridIPE) then
        CALL check2GridIPE(unitCheckIPE, NPTS, NMP, &
                           glat_plasma_3d, glond_plasma_3d, pz_plasma_3d) 

        CLOSE (unitCheckIPE)

    endif ! (debugGridIPE)


    !--------------------------------------
    ! Define name of the geo to mag file
    !--------------------------------------
    geoToMagFileName = TRIM(staticFileDir)//'/Geographic_to_Magnetic_91_20.2000.0.format'


    !------------------------------------------------
    ! Open & read geoToMagFileName
    ! This subroutine is in moduleThermosphere :
    !------------------------------------------------
    call calculate_magnetic_parameters_using_apex( geoToMagFileName, &
               B_magnitude_apex_nT, B_dip_angle_apex_degrees, &
               Magnetic_latitude_degrees, Magnetic_longitude_degrees)

	      
    !--------------------------------------------------------------------
    ! Calculate time for initializing the thermosphere
    !--------------------------------------------------------------------
    Universal_Time_seconds = start_time 

    IF ( .NOT. switches%electro ) THEN 
        OPEN (41, FILE=TRIM(staticFileDir)//&
                            &'jicamarca_zonal_drifts', STATUS='OLD')

        !-----------------------------------------------
        ! Inputs are f107, nday, the rest are outputs
        ! This subroutine reads from file unit #41
        !-----------------------------------------------
        call low_lat_efield (exns, eyns, ezns, plvu, zonal, f107, nday)
    ENDIF 

    !---------------------------------
    ! Initialize nnloop and idump_gt
    !---------------------------------
    nnloop = 0
    idump_gt = 0

    !---------------------------------------------------------
    ! Initialize counter for smoothing and for composition
    !---------------------------------------------------------
    nn_smoothing_counter = 0    !lrm20120323
    nn_composition_counter = 0  !lrm20120323


    CALL GT_thermosphere_INIT( &
                 GT_input_dataset, &  ! startup file                               ! INPUT
                 GT_output_dataset, & ! output file                                ! INPUT
                 nday, &                                                           ! INPUT
                 Universal_Time_seconds, &                                         ! INPUT
                 solar_declination_angle_radians, &                                ! OUTPUT
                 hough11 , hough22 , hough23 , hough24 , hough25, &                ! OUTPUT
                 amplitude, &                                                      ! INPUT
                 tidalPhase, &
                 wind_southwards_ms1_FROM_GT, &  ! read from startup file             ! OUTPUT
                 wind_eastwards_ms1_FROM_GT, & ! read from startup file               ! OUTPUT
                 wvz_FROM_GT, &  ! upward nuetral wind                                ! OUTPUT
                 rmt_FROM_GT, & ! composition                                         ! OUTPUT
                 Temperature_K_FROM_GT, & ! nuetral temperature in 3D pressure coords ! OUTPUT
                 ht_FROM_GT, & ! height at each grid point in 3D                        ! OUTPUT
                 Ne_density_FOR_GT, &    !ionospheric electron density interpolated to pressure heights OUTPUT
                 Oplus_density_FOR_GT, & !ionospheric O+ density interpolated to pressure heights OUTPUT
                 NOplus_density_FOR_GT, &!ionospheric NO+ density interpolated to pressure heights OUTPUT
                 O2plus_density_FOR_GT, &!ionospheric O2+ density interpolated to pressure heights OUTPUT
                 Ti_Oplus_FOR_GT )       !ionospheric Ti O+ temperature interpolated to pressure heights OUTPUT


      print *,'driver_ipe_gt.3d : ****************************************************************'
      print *,'driver_ipe_gt.3d : USING COLD START VALUES FOR IONOSPHERE *************************'
      print *,'driver_ipe_gt.3d : ****************************************************************'

      Oplus_density_FOR_GT(:,:,:) = 1.e09
      NOplus_density_FOR_GT(:,:,:) = 1.e03
      O2plus_density_FOR_GT(:,:,:) = 1.e03
      Ti_Oplus_FOR_GT(:,:,:) = 1000.
      Ne_density_FOR_GT(:,:,:) = 2.e09

      !--------------------------------------------------
      ! Initialize all of the GIP switches to be .FALSE.
      !--------------------------------------------------
      GIP_switches(:) = .FALSE.

      !--------------------------------------------------
      ! then set the the ones that are used..
      ! THESE ARE SET IN THE NAME LIST NOW ********
      !--------------------------------------------------

      GIP_switches(5) = switches%External_model_provides_NO_N4S_densities
      GIP_switches(6) = switches%External_model_provides_low_lat_E_fields
      GIP_switches(7) = switches%input_Auroral_production_is_single_overall_rate


      !-----------------------------------------------------------
      ! Set up lat, longs, used for interpolation to the IPE grid
      !-----------------------------------------------------------
      do ii = 1 , GT_lon_dim
         therm_model_geo_long_deg(ii) = (ii - 1) * 18.
      enddo

      do ii = 1 , GT_lat_dim
         therm_model_geo_lat_deg(ii) = (ii - 46) * 2.
      enddo


      !-------------------------------------------------
      ! Debug the results from the GT_thermosphere_INIT
      !-------------------------------------------------
      If (debugThermo) then

	   print *,'driver_ipe_gt.3d : AFTER CALLING GT_thermosphere_INIT.............'

           CALL checkThermo(debugDir, debugThermoFileName, unitCheckThermo, &
                            GT_input_dataset, GT_output_dataset, nday, &
                            Universal_Time_seconds, solar_declination_angle_radians,  &
                            houghSize, hough11, hough22, hough23, hough24, hough25, &
                            amplitude, tidalPhase, GT_ht_dim, GT_lat_dim, GT_lon_dim, &
                            wind_southwards_ms1_FROM_GT, &
                            wind_eastwards_ms1_FROM_GT, wvz_FROM_GT, rmt_FROM_GT, &
                            Temperature_K_FROM_GT, ht_FROM_GT, Ne_density_FOR_GT)

           STOP


       endif ! debugThermo

     !----------------------------------------------
     ! Read in the fixed grid ionospheric grid
     ! (only need to do this once)
     !----------------------------------------------
     CALL readIPEtoGeoGrid( giptogeoFileName)


endif ! sw_neutral == 'GT' 


!------------------------------------------------
! Open IPE species files
! One for each species, with all time steps in it
!------------------------------------------------




!-----------------------------------------------
!  Ionospheric Loop :  time_loop is in seconds
!-----------------------------------------------
time_loop: DO utime = (start_time + GT_timestep_in_seconds), stop_time, time_step
  !print *, 'driver_ipe_gt.3d : utime  =  ', utime


  !**************************************************************************************
  ! somewhere in the loop read startup00.ascii (Oplus)  and all the other species
  ! in a separate loop  (want to read all the species from IPE at a certain time)
  ! To read binary files, just remove the format statement
  !**************************************************************************************
  ! read the O plus data (ascii file) from startup file
  !O+ density
    

  !------------------------------
  ! Read the IPE start up files
  !------------------------------
  READ (UNIT=startUpFiles(1)%unit,FMT="(20E12.4)" ) Oplus_density_from_IPE
  READ (UNIT=startUpFiles(2)%unit,FMT="(20E12.4)" ) Hplus_density_from_IPE
  READ (UNIT=startUpFiles(3)%unit,FMT="(20E12.4)" ) Heplus_density_from_IPE
  READ (UNIT=startUpFiles(4)%unit,FMT="(20E12.4)" ) Nplus_density_from_IPE
  READ (UNIT=startUpFiles(5)%unit,FMT="(20E12.4)" ) NOplus_density_from_IPE
  READ (UNIT=startUpFiles(6)%unit,FMT="(20E12.4)" ) O2plus_density_from_IPE
  READ (UNIT=startUpFiles(7)%unit,FMT="(20E12.4)" ) N2plus_density_from_IPE
  READ (UNIT=startUpFiles(8)%unit,FMT="(20E12.4)" ) Oplus2D_density_from_IPE
  READ (UNIT=startUpFiles(9)%unit,FMT="(20E12.4)" ) Oplus2P_density_from_IPE
  READ (UNIT=startUpFiles(10)%unit,FMT="(20E12.4)" ) Te_from_IPE  
  READ (UNIT=startUpFiles(11)%unit,FMT="(20E12.4)" ) Ti_from_IPE  

! Not using unit 12 right now lrm20120831

! WILL NEED TO READ VELOCITY FOR OPLUS AND HPLUS AND INTERPOLATE IN THE FUTURE may 1, 2012

  print *,'driver_ipe_gt.3d '
  print *,'driver_ipe_gt.3d : after reading all ipe startup files.......'

  !--------------------------------------------------------------------------
  !Ne_density_from_IPE  = sum of all above ion densities  (not temperatures)
  !--------------------------------------------------------------------------
  Ne_density_from_IPE  = Oplus_density_from_IPE + Hplus_density_from_IPE + &
                         Heplus_density_from_IPE + Nplus_density_from_IPE + &
                         NOplus_density_from_IPE + O2plus_density_from_IPE + &
                         N2plus_density_from_IPE + Oplus2D_density_from_IPE + &
                         Oplus2P_density_from_IPE
  
  If (debugStartupIPE) then 

      print *,'driver_ipe_gt.3d : Oplus_density_from_IPE(1:5,1) = ', Oplus_density_from_IPE(1:5,1)
      print *,'driver_ipe_gt.3d : Hplus_density_from_IPE(1:5,1) = ', Hplus_density_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : Heplus_density_from_IPE(1:5,1) = ', Heplus_density_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : Nplus_density_from_IPE(1:5,1) = ', Nplus_density_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : NOplus_density_from_IPE(1:5,1) = ', NOplus_density_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : O2plus_density_from_IPE(1:5,1) = ', O2plus_density_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : N2plus_density_from_IPE(1:5,1) = ', N2plus_density_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : Oplus2D_density_from_IPE(1:5,1) = ', Oplus2D_density_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : Oplus2P_density_from_IPE(1:5,1) = ', Oplus2P_density_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : Te_from_IPE(1:5,1) = ', Te_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : Ti_from_IPE(1:5,1) = ', Ti_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d : Ne_density_from_IPE(1:5,1) = ', Ne_density_from_IPE(1:5,1)
       print *,'driver_ipe_gt.3d '

       print *,'driver_ipe_gt.3d : inGIP1d(1:NLP) = ', inGIP1d(1:NLP)
       print *,'driver_ipe_gt.3d : isGIP1d(1:NLP) = ', isGIP1d(1:NLP)


     ! Note :
     !JMIN_IN(lp) <=== inGIP1d(1:NLP)
     !JMAX_IS(lp) <=== isGIP1d(1:NLP)


     ! For checking ...........
     !do mm = 1, NMP
!
     !   do ll = 1,6
     !       jj = inGIP1d(ll)  
     !       kk = isGIP1d(ll)  

     !       print *, jj, kk

     !       do ii = kk, kk+3
     !          print *, ii, ll, mm, Oplus_density_from_IPE(ii,mm)
     !       end do !ii
     !   end do !ll

     !end do !mm



     print *,'driver_ipe_gt.3d : MINVAL(Oplus_density_from_IPE) = ', MINVAL(Oplus_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(Oplus_density_from_IPE) = ', MAXVAL(Oplus_density_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(Hplus_density_from_IPE) = ', MINVAL(Hplus_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(Hplus_density_from_IPE) = ', MAXVAL(Hplus_density_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(Heplus_density_from_IPE) = ', MINVAL(Heplus_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(Heplus_density_from_IPE) = ', MAXVAL(Heplus_density_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(Nplus_density_from_IPE) = ', MINVAL(Nplus_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(Nplus_density_from_IPE) = ', MAXVAL(Nplus_density_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(NOplus_density_from_IPE) = ', MINVAL(NOplus_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(NOplus_density_from_IPE) = ', MAXVAL(NOplus_density_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(O2plus_density_from_IPE) = ', MINVAL(O2plus_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(O2plus_density_from_IPE) = ', MAXVAL(O2plus_density_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(N2plus_density_from_IPE) = ', MINVAL(N2plus_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(N2plus_density_from_IPE) = ', MAXVAL(N2plus_density_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(Oplus2D_density_from_IPE) = ', MINVAL(Oplus2D_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(Oplus2D_density_from_IPE) = ', MAXVAL(Oplus2D_density_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(Oplus2P_density_from_IPE) = ', MINVAL(Oplus2P_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(Oplus2P_density_from_IPE) = ', MAXVAL(Oplus2P_density_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(Te_from_IPE) = ', MINVAL(Te_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(Te_from_IPE) = ', MAXVAL(Te_from_IPE)

     print *,'driver_ipe_gt.3d : MINVAL(Ti_from_IPE) = ', MINVAL(Ti_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(Ti_from_IPE) = ', MAXVAL(Ti_from_IPE)


     print *,'driver_ipe_gt.3d : MINVAL(Ne_density_from_IPE) = ', MINVAL(Ne_density_from_IPE)
     print *,'driver_ipe_gt.3d : MAXVAL(Ne_density_from_IPE) = ', MAXVAL(Ne_density_from_IPE)


     STOP

  endif ! debugStartupIPE

  If (debugIonoInterp) then 

      !20120919lrm

      write(unitCheckIonoInterpBefore,*) Oplus_density_from_IPE
      write(unitCheckIonoInterpBefore + 1,*) Hplus_density_from_IPE
      write(unitCheckIonoInterpBefore + 2,*) Nplus_density_from_IPE
      write(unitCheckIonoInterpBefore + 3,*) NOplus_density_from_IPE
      write(unitCheckIonoInterpBefore + 4,*) O2plus_density_from_IPE
      write(unitCheckIonoInterpBefore + 5,*) N2plus_density_from_IPE
      write(unitCheckIonoInterpBefore + 6,*) Ne_density_from_IPE
      write(unitCheckIonoInterpBefore + 7,*) Te_from_IPE
      write(unitCheckIonoInterpBefore + 8,*) Ti_from_IPE

  end if

  !---------------------------------------------------------------
  ! Check the cpu time to get timing of the interface subroutine
  !---------------------------------------------------------------
  CALL CPU_TIME(startTime)


  !-------------------------------------------------------
  ! 1st interpolation :  Interpolate from flux tube grid
  ! to the fixed height grid for all IPE species
  !-------------------------------------------------------
  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   Oplus_density_from_IPE, &      ! inputs
                                   Oplus_high_res_fixed)   ! output

  CALL CPU_TIME(endTime)

  !------------------------------------------------------------------------
  ! Add up how much time we've been spent in this routine, then avg at end
  !------------------------------------------------------------------------
  timeFluxGridtoFixedGrid = (endTime-startTime) + timeFluxGridtoFixedGrid
  WRITE(*,*) "cpu_time for INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO, Oplus   : ",(endTime-startTime)




  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   Hplus_density_from_IPE, &      ! inputs
                                   Hplus_high_res_fixed)   ! output


  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   Nplus_density_from_IPE, &      ! inputs
                                   Nplus_high_res_fixed)   ! output



  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   NOplus_density_from_IPE, &      ! inputs
                                   NOplus_high_res_fixed)   ! output



  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   O2plus_density_from_IPE, &      ! inputs
                                   O2plus_high_res_fixed)   ! output



  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   N2plus_density_from_IPE, &      ! inputs
                                   N2plus_high_res_fixed)   ! output


				   
  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   Ne_density_from_IPE, &      ! inputs
                                   Ne_high_res_fixed)   ! output



  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   Te_from_IPE, &      ! inputs
                                   Te_high_res_fixed)   ! output



  CALL INTERFACE__MID_LAT_IONOSPHERE_to_FIXED_GEO( &
                                   NPTS, NMP, NLP, &     ! inputs
                                   nFixedGridIonoHeights, nFixedGridIonoLats, nFixedGridIonoLons, & ! inputs
                                   Ti_from_IPE, &      ! inputs
                                   Ti_high_res_fixed)   ! output		   
  If (debugIonoInterp) then 
      !print *,'driver_ipe_gt.3d : MINVAL(Oplus_high_res_fixed) = ', MINVAL(Oplus_high_res_fixed)
      !nm20120607dbg
      !write(9999,*) Oplus_high_res_fixed   ! output
      write(unitCheckIonoInterpAfter,*) Oplus_high_res_fixed
      write(unitCheckIonoInterpAfter+1,*) Hplus_high_res_fixed
      write(unitCheckIonoInterpAfter+2,*) Nplus_high_res_fixed
      write(unitCheckIonoInterpAfter+3,*) NOplus_high_res_fixed
      write(unitCheckIonoInterpAfter+4,*) O2plus_high_res_fixed
      write(unitCheckIonoInterpAfter+5,*) N2plus_high_res_fixed
      write(unitCheckIonoInterpAfter+6,*) Ne_high_res_fixed
      write(unitCheckIonoInterpAfter+7,*) Te_high_res_fixed
      write(unitCheckIonoInterpAfter+8,*) Ti_high_res_fixed
  end if


  IF ( sw_neutral == 'GT' ) then

   !------------------------------------------
   ! Loop here for the thermosphere time step
   !------------------------------------------
   thermosphereLoop : DO gtLoopTime = utime, utime + (time_step-GT_timestep_in_seconds), &
                                      GT_timestep_in_seconds !------------------------

      nnloop = MOD(nnloop + 1, number_of_GT_time_steps_in_24_hours)
      IF ( nnloop .EQ. 0 ) nnloop = number_of_GT_time_steps_in_24_hours

      Universal_Time_seconds = gtLoopTime

      nn = MOD(nnloop, number_of_GT_time_steps_in_24_hours)
      IF ( nn .EQ. 0 ) nn = number_of_GT_time_steps_in_24_hours

      print *, 'driver_ipe_gt.3d : utime, Universal_Time_seconds, nn, nnloop  =  ', &
                                   utime, Universal_Time_seconds, nn, nnloop

      !---------------------
      ! increment counters
      !---------------------
      nn_smoothing_counter = nn_smoothing_counter + 1
      nn_composition_counter = nn_composition_counter + 1

      !----------------------------------------------------
      ! Set Foster parameters
      ! These will change in time in future versions
      !----------------------------------------------------
      Foster_level = 5

      Foster_power = 130.  ! Note: this is only used if Foster_level (above) is 10.

      Foster_Efield_amplification = 1.3

      CALL FOSTER(exns, eyns, Foster_level, Foster_power, &
                  Foster_Efield_amplification, Foster_potential)

      !------------------------------------------------------------------
      ! The interface from fixed grid (fixed height) to GT pressure GRID
      !-------------------------------------------------------------------
      ! from res-thermo-02:/ipe/lmayer/gt_gip_v1.2/GIP_ionosphere_plasmasphere.f90 -
      !  - was INTERFACE__GIP_to_thermosphere
      ! Created an independent subroutine out of the module - 20120501lrm
      ! - renamed this ipe fixed grid to thermosphere - 20120501lrm
      ! call INTERFACE__FIXED_GRID_to_THERMO ( &
      !    thermospheric_model_name , therm_model_ht_dim , therm_model_lat_dim , therm_model_lon_dim, &
      !    ne_high_res_fixed, Oplus_high_res_fixed, hplus_high_res_fixed, &
      !    noplus_high_res_fixed, o2plus_high_res_fixed, &
      !    n2plus_high_res_fixed, nplus_high_res_fixed, &
      !    Te_high_res_fixed, Ti1_high_res_fixed, Ti2_high_res_fixed, &
      !    therm_model_geo_long, therm_model_geo_lat, therm_model_ht_m, &
      !    therm_model_Ne_density, therm_model_oplus_density, therm_model_hplus_density, &
      !    therm_model_noplus_density, therm_model_o2plus_density, &
      !    therm_model_n2plus_density, therm_model_nplus_density, &
      !    therm_model_Te, therm_model_Ti1, therm_model_Ti2)

      ! from res-thermo-02 :
      !   call INTERFACE__GIP_to_thermosphere ( &
      !         thermospheric_model_name , therm_model_ht_dim , therm_model_lat_dim , therm_model_lon_dim, &
      !         ne_high_res_fixed,oplus_high_res_fixed,hplus_high_res_fixed, &
      !         noplus_high_res_fixed,o2plus_high_res_fixed, &
      !         n2plus_high_res_fixed,nplus_high_res_fixed, &
      !         Te_high_res_fixed,Ti1_high_res_fixed,Ti2_high_res_fixed, &
      !         therm_model_geo_long, therm_model_geo_lat, therm_model_ht_m, &
      !         therm_model_Ne_density,therm_model_oplus_density,therm_model_hplus_density, &
      !         therm_model_noplus_density,therm_model_o2plus_density, &
      !         therm_model_n2plus_density,therm_model_nplus_density, &
      !         therm_model_Te,therm_model_Ti1,therm_model_Ti2)

      !----------------
      ! Placeholders :
      !----------------
      !Hplus_high_res_fixed = -999999.  ! This is not used yet ****  lrm20120509
      !N2plus_high_res_fixed = -999999.  ! This is not used yet ****  lrm20120509
      !Nplus_high_res_fixed = -999999.  ! This is not used yet ****  lrm20120509
      !Te_high_res_fixed = -999999.  ! This is not used yet ****  lrm20120509

      !print *,'driver : Before INTERFACE__FIXED_GRID_to_THERMO '
      !print *,'driver : Ne_high_res_fixed(:,10,10) = ', Ne_high_res_fixed(:,10,10)

      !------------------------------------------------------
      ! Time the fixed grid to pressure grid interpolation
      !------------------------------------------------------
      CALL cpu_time(startTime)

      ! This is in the loop b/c of changing heights of the pressure grid
      CALL INTERFACE__FIXED_GRID_to_THERMO ( &
         thermospheric_model_name , GT_ht_dim , GT_lat_dim , GT_lon_dim , &  ! inputs
         Ne_high_res_fixed, Oplus_high_res_fixed, Hplus_high_res_fixed, &    ! inputs
         NOplus_high_res_fixed, O2plus_high_res_fixed, &                     ! inputs
         N2plus_high_res_fixed, Nplus_high_res_fixed, &                      ! inputs
         Te_high_res_fixed, Ti_high_res_fixed, Ti_high_res_fixed, &          ! inputs
         therm_model_geo_long_deg, therm_model_geo_lat_deg, &                ! inputs
         !Altitude_m_FOR_IPE, &                                               ! inputs   *** this is not set yet***
         Ht_FROM_GT, &                                               ! inputs   ht, lat, lon
         Ne_density_FOR_GT, Oplus_density_FOR_GT, Hplus_density_FOR_GT, &    ! outputs
         NOplus_density_FOR_GT, O2plus_density_FOR_GT, &                     ! outputs
         N2plus_density_FOR_GT, Nplus_density_FOR_GT, &                      ! outputs
         Te_FOR_GT, Ti_Oplus_FOR_GT, Ti_Oplus_FOR_GT)                        ! outputs

      CALL cpu_time(endTime)

      !------------------------------------------------------------------------
      ! Add up how much time we've been spent in this routine, then avg at end
      !------------------------------------------------------------------------
      timeFixedGridtoPressureGrid = (endTime-startTime) + timeFixedGridtoPressureGrid
      WRITE(*,*) "cpu_time for  INTERFACE__FIXED_GRID_to_THERMO  : ",(endTime-startTime)


         !--------------------------------------------------------
         ! Write out the interpolated values on the pressure grid
         !--------------------------------------------------------
         if ((debugIonoFixedtoPressure) .and. (idump_gt == 1)) then

            write(unitCheckPressureInterp,*) Oplus_density_for_GT
            write(unitCheckPressureInterp + 1,*) Hplus_density_for_GT
            write(unitCheckPressureInterp + 2,*) Nplus_density_for_GT
            write(unitCheckPressureInterp + 3,*) NOplus_density_for_GT
            write(unitCheckPressureInterp + 4,*) O2plus_density_for_GT
            write(unitCheckPressureInterp + 5,*) N2plus_density_for_GT
            write(unitCheckPressureInterp + 6,*) Ne_density_for_GT
            write(unitCheckPressureInterp + 7,*) Te_for_GT
            write(unitCheckPressureInterp + 8,*) Ti_Oplus_for_GT

            ! write out the pressure grid heights
            write(unitPressureGridHeight,*) ht_FROM_GT

           !--------------------------------------------------------------------
           ! Calculate and write out the average height for each pressure level
           !--------------------------------------------------------------------
           write(unitCheckPressureHeight, FMT="(I12)" ) gtLoopTime
           do ii = 1, GT_ht_dim

              !avgHeightofPressure = Ht_FROM_GT(ii,:,:)/(GT_lat_dim*GT_lon_dim)

              write(unitCheckPressureHeight, FMT="(E12.4)") SUM(Ht_FROM_GT(ii,:,:))/(GT_lat_dim*GT_lon_dim)
              !write(unitCheckPressureHeight, FMT="(A1)") " "

           enddo

        end if


        !----------------
        ! Placeholders :
        !----------------
        !Hplus_density_FOR_GT = -999999.  ! This is not used yet ****  20120508lrm
        !N2plus_density_FOR_GT = -999999.  ! This is not used yet ****  20120508lrm
        !Nplus_density_FOR_GT = -999999.  ! This is not used yet ****  20120508lrm
        !Te_FOR_GT = -999999.  ! This is not used yet ****  20120508lrm


        !--------------------------------------------
        ! Variables for GT - 20120501lrm
        !--------------------------------------------
        !therm_geo_long_input : therm_model_geo_long_deg
        !therm_geo_lat : therm_model_geo_lat_deg
        !therm_Z : ht_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim) : Altitude_m_FOR_IPE(GT_ht_dim, GT_lon_dim, GT_lat_dim)

        !therm_Ne_density : Ne_density_FOR_GT
        !therm_oplus_density : Oplus_density_FOR_GT
        !therm_hplus_density ????   keep  - not sure if we want to use this ???
        !therm_noplus_density : NOplus_density_FOR_GT
        !therm_o2plus_density : O2plus_density_FOR_GT
        !therm_n2plus_density : ???? keep
        !therm_nplus_density : ???? keep
        !therm_Te : ???? keep  - should be used later
        !therm_Ti1 : Ti_Oplus_FOR_GT - yes, use Ti
        !therm_Ti2 : only 1 Ti from IPE

        !print *,'driver : After INTERFACE__FIXED_GRID_to_THERMO '
        !print *,'driver : Ne_density_FOR_GT(:,1,1) = ',Ne_density_FOR_GT(:,1,1)


        ! will need to add heating rates from IPE - Tim will add heating rates to GT

        !----------------------------------------
        ! Dump out the thermosphere values ???
        !----------------------------------------
        if (nnloop .eq. nnstop) idump_gt = 1

        !---------------------------------------------
        ! Now call the thermosphere ................
        !---------------------------------------------

        ! This will reset the nn_smoothing_counter if it equals the smoothing frequency
        call GT_thermosphere( &
                      GT_output_dataset, &
                      GT_timestep_in_seconds, &
                      idump_gt, &
                      solar_declination_angle_radians, &
                      nn, nnloop, &
                      Universal_Time_seconds, &
                      nn_smoothing_counter, &
                      smoothingFrequency, &
                      nn_composition_counter, &
                      neutralCompositionFrequency, &
                      hough11, hough22, hough23, hough24, hough25, &
                      amplitude, &
                      tidalPhase, &
                      Ne_density_FOR_GT, &      ! input/output
                      Oplus_density_FOR_GT, &   ! input
                      NOplus_density_FOR_GT, &  ! input
                      O2plus_density_FOR_GT, &  ! input
                      Ti_Oplus_FOR_GT, &        ! input
                      exns, eyns, ezns, &
                      B_dip_angle_apex_degrees, B_magnitude_apex_nT, &
                      Magnetic_latitude_degrees, Magnetic_longitude_degrees, &
                      Foster_level, &  ! newl
                      Foster_power, &  ! gw
                      f107, &
                      emaps, cmaps, profile, &
                      wind_southwards_ms1_FROM_GT, &  ! output
                      wind_eastwards_ms1_FROM_GT, &  ! output
                      wvz_FROM_GT, &                  ! output
                      rmt_FROM_GT, &  ! output
                      Temperature_K_FROM_GT, &  ! output
                      ht_FROM_GT, &  ! output
                      O_density_FROM_GT, &  ! output
                      O2_density_FROM_GT, &  ! output
                      N2_density_FROM_GT, &  ! output
                      qion3d)  ! output


       !-------------------------------------------------------------------
       ! If in Thermospheric debug mode, write out values to a file
       ! & do some prints to the screen
       !-------------------------------------------------------------------
       if ((debugThermo).and. (idump_gt == 1)) then

          WRITE (unitCheckThermo, '(A)' ) ' '
          WRITE(unitCheckThermo,'(A)') 'driver_ipe_gt.3d : '
          WRITE(unitCheckThermo,'(A, i5)') 'driver_ipe_gt.3d : utime =  ' , utime
          WRITE(unitCheckThermo,'(A, i5)') 'driver_ipe_gt.3d : nn, nnloop=  ' , nn, nnloop
          write(unitCheckThermo,'(A, 1F20.7)') 'driver_ipe_gt.3d : Universal_Time_seconds = ' , &
                                                Universal_Time_seconds
          write(unitCheckThermo,'(A)') 'driver_ipe_gt.3d : ----------------------------------------------------'

          CALL checkThermoArrays(unitCheckThermo, &
                                 GT_ht_dim, GT_lat_dim, GT_lon_dim, &
                                 wind_southwards_ms1_FROM_GT, &
                                 wind_eastwards_ms1_FROM_GT, wvz_FROM_GT, rmt_FROM_GT, &
                                 Temperature_K_FROM_GT, ht_FROM_GT)

        endif ! debugThermo


        littleLoop = littleLoop + 1

      END DO thermosphereLoop ! --------------------------------------------------------------------------------------------------



        if (debugThermoInterp) then

           !print *,'driver : setting GT densities to constants for testing ********************' 
           !print *,'driver : setting GT densities to constants for testing ********************' 
           !print *,'driver : setting GT densities to constants for testing ********************' 


           !do ii = 1, gt_ht_dim
           !   O_density_FROM_GT(ii,:,:) = ii
           !end do

           !O_density_FROM_GT = iitime

           !iitime = iitime + 1


           !O_density_FROM_GT = 1.0
           !O2_density_FROM_GT = 2.0
           !N2_density_FROM_GT = 3.0

           ! write out arrays to file for plotting
           write(unitCheckThermoInterpBefore,*) wind_southwards_ms1_FROM_GT
           write(unitCheckThermoInterpBefore + 1,*) wind_eastwards_ms1_FROM_GT
           write(unitCheckThermoInterpBefore + 2,*) wvz_FROM_GT
           write(unitCheckThermoInterpBefore + 3,*) rmt_FROM_GT
           write(unitCheckThermoInterpBefore + 4,*) Temperature_K_FROM_GT
           write(unitCheckThermoInterpBefore + 5,*) O_density_FROM_GT
           write(unitCheckThermoInterpBefore + 6,*) O2_density_FROM_GT
           write(unitCheckThermoInterpBefore + 7,*) N2_density_FROM_GT
           write(unitCheckThermoInterpBefore + 8,*) qion3d

          !--------------------------------------------------------------------
          ! Calculate and write out the average height for each pressure level
          !--------------------------------------------------------------------
          write(unitCheckPressureHeightThermo, FMT="(I12)" ) gtLoopTime
         do ii = 1, GT_ht_dim

            write(unitCheckPressureHeightThermo, FMT="(E12.4)") SUM(Ht_FROM_GT(ii,:,:))/(GT_lat_dim*GT_lon_dim)

         enddo

      endif ! debugThermoInterp






      !-------------------------------------------------------------------------
      ! For interpolation to IPE grid, switch (ht, lat, lon) to (ht, lon, lat)
      !-------------------------------------------------------------------------
      do ii = 1 , GT_ht_dim
         do jj = 1 , GT_lon_dim
            do kk = 1 , GT_lat_dim

               Altitude_m_FOR_IPE(ii,jj,kk) = ht_FROM_GT(ii,kk,jj)
               Vn_Southwards_ms1_FOR_IPE(ii,jj,kk) = wind_southwards_ms1_FROM_GT(ii,kk,jj)
               Vn_Eastwards_ms1_FOR_IPE(ii,jj,kk) = wind_eastwards_ms1_FROM_GT(ii,kk,jj)
               Vn_Upwards_ms1_FOR_IPE(ii,jj,kk) = wvz_FROM_GT(ii,kk,jj)
               Tn_K_FOR_IPE(ii,jj,kk) = temperature_K_FROM_GT(ii,kk,jj)
               O_density_m3_FOR_IPE(ii,jj,kk) = O_density_FROM_GT(ii,kk,jj)
               O2_density_m3_FOR_IPE(ii,jj,kk) = O2_density_FROM_GT(ii,kk,jj)
               N2_density_m3_FOR_IPE(ii,jj,kk) = N2_density_FROM_GT(ii,kk,jj)
               elx_FOR_IPE(ii,jj,kk) = 0.  ! = elx(m,l) ****** elx IS NOT DEFINED IN TUCAN_TIME ***********
               ely_FOR_IPE(ii,jj,kk) = 0.  ! = ely(m,l) ****** ely IS NOT DEFINED IN TUCAN_TIME ***********
               qion3d_FOR_IPE(ii,jj,kk) = qion3d(ii,kk,jj)

            enddo
         enddo
      enddo


   
      !***************************************************
      !!! HN_m3,HE_m3 still need to be obtained from MSIS
      ! call MSIS to get these values
      !***************************************************
      ! I could use the nday from the INP.inp file ****************
      !iday = nday
      !iyd = 99000 + iday
      !call gtd7(iyd,sec,alt,glat,glon,stl,f107a_msis,f107d_msis,ap_msis,mass,d,t)

      ! convert from cm-3 to m-3
      !he_density_m3 (i) = d(1)/M3_TO_CM3 !*1.e6
      !h_density_m3  (i) = d(7)/M3_TO_CM3
      !n4s_density_m3(i) = d(8)/M3_TO_CM3

      !-------------------------------------------------------------------------
      ! time how long it takes to do pressure grid to fixed grid interpolation
      !-------------------------------------------------------------------------
      CALL cpu_time(startTime)

      !------------------------------------------------------
      ! interpolate from pressure grid to fixed height grid
      !------------------------------------------------------
        call INTERFACE__thermosphere_to_FIXED_GEO ( &
           GIP_switches, &
           thermospheric_model_name , &
           GT_ht_dim, &                    ! was : therm_model_ht_dim
           GT_lat_dim, &                   ! was : therm_model_lat_dim
           GT_lon_dim, &                   ! was : therm_model_lon_dim
           O_density_m3_FOR_IPE, & !  O_density_FROM_GT,  &             ! was : therm_model_o_density
           O2_density_m3_FOR_IPE, & ! O2_density_FROM_GT, &             ! was : therm_model_o2_density
           N2_density_m3_FOR_IPE, & ! N2_density_FROM_GT, &             ! was : therm_model_n2_density
           !! therm_model_NO_density, &      ! not needed now
           !! therm_model_N4S_density, &     ! not needed now
           !! therm_model_N2D_density, &     ! not needed now
           Tn_K_FOR_IPE, & ! Temperature_K_FROM_GT, &        ! was : therm_model_Tn
           Vn_Eastwards_ms1_FOR_IPE, & ! wind_eastwards_ms1_FROM_GT, &  ! was : therm_model_Vy
           Vn_Southwards_ms1_FOR_IPE, &  ! wind_southwards_ms1_FROM_GT, &   ! was : therm_model_Vx
           Vn_Upwards_ms1_FOR_IPE, &  ! wind_upwards_ms1_from_gt,  added lrm20121108
           qion3d_FOR_IPE, & ! was : therm_model_qion3d    ! input
           elx_FOR_IPE, & ! exns, &  ****                       ! was : therm_model_elx
           ely_FOR_IPE, & ! eyns, &  ****                      ! was : therm_model_ely
           
           !! AURORA VARIABLES NOT NEEDED FOR IPE YET lrm20110929
           !! therm_model_qo2p_aurora, therm_model_qop_aurora, therm_model_qn2p_aurora, &
           !! therm_model_qnp_aurora, therm_model_qtef_aurora, &

           therm_model_geo_long_deg, &  ! was : therm_model_geo_long
           therm_model_geo_lat_deg, &   ! was : therm_model_geo_lat 
           Altitude_m_FOR_IPE, &        ! was : therm_model_ht_m, &

           O_density_fixed_ht, O2_density_fixed_ht, N2_density_fixed_ht, &     ! output


           !! NO_density_fixed_ht, &    ! output, not needed now
           !! N4S_density_fixed_ht, &   ! output, not needed now
           !! N2D_density_fixed_ht, &   ! output, not needed now


           V_East_FixedHeight, V_South_FixedHeight, V_Upward_FixedHeight, tts_fixed_ht, &             ! output
           qion3d_fixed_ht, elx_fixed_ht, ely_fixed_ht)                        ! output

           !! AURORA VARIABLES NOT NEEDED FOR IPE YET lrm20110929
           !! qo2p_aurora_fixed_ht, qop_aurora_fixed_ht, qn2p_aurora_fixed_ht, &  ! output
           !! qnp_aurora_fixed_ht, qtef_aurora_fixed_ht)                          ! output


           CALL cpu_time(endTime)
           timePressureGridtoFixedGrid = (endTime - startTime) + timePressureGridtoFixedGrid 
           WRITE(*,*) "cpu_time for INTERFACE__thermosphere_to_FIXED_GEO   : ",(endTime-startTime)

           !-------------------------------------------------------------------------
           ! Write out results of the interpolation to an ascii file for examination
           !-------------------------------------------------------------------------
           if (debugThermoInterp) then

              print *,'driver_ipe_gt.3d : Writing out the interpolated values to a file..........'

              CALL checkInterp(debugDir, debugThermoInterpFileName, unitCheckThermoInterp, &
                       nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim, &
                       O_density_fixed_ht, O2_density_fixed_ht, N2_density_fixed_ht, &
                       V_East_FixedHeight, V_South_FixedHeight, V_Upward_FixedHeight, &
                       tts_fixed_ht, qion3d_fixed_ht, elx_fixed_ht, &
                       ely_fixed_ht)

               !print *,'driver_ipe_gt.3d : STOPPING.............'
               !STOP
           endif

           !------------------------------------------------------
           ! Write out pressure to fixed grid results
           !------------------------------------------
  If (debugThermoInterp) then 
      write(unitCheckThermoInterpAfter,*) V_South_FixedHeight
      write(unitCheckThermoInterpAfter+1,*) V_East_FixedHeight
      write(unitCheckThermoInterpAfter+2,*) V_Upward_FixedHeight
      write(unitCheckThermoInterpAfter+3,*) ! should be rmt , but the composition is not being used ***
      write(unitCheckThermoInterpAfter+4,*) tts_fixed_ht
      write(unitCheckThermoInterpAfter+5,*) O_density_fixed_ht
      write(unitCheckThermoInterpAfter+6,*) O2_density_fixed_ht
      write(unitCheckThermoInterpAfter+7,*) N2_density_fixed_ht
      write(unitCheckThermoInterpAfter+8,*) qion3d_fixed_ht
  end if

  !---------------------------------------------------------------------
  ! Time how long it takes to interpolate fixed grid to flux tube grid
  !---------------------------------------------------------------------
  CALL cpu_time(startTime)


  !--------------------------------------------
  ! Convert fixed height grid to flux tube grid
  !--------------------------------------------


  call INTERFACE__FIXED_GEO_to_IONOSPHERE( &
        therm_model_geo_long_deg, & ! was : geo_grid_longitudes_degrees, &
        therm_model_geo_lat_deg, &  ! was : geo_grid_latitudes_degrees, &
        O_density_fixed_ht, &
        O2_density_fixed_ht, &
        N2_density_fixed_ht, &

        !! NO_density_fixed_ht, &    not needed now
        !! N4S_density_fixed_ht, &   not needed now
        !! N2D_density_fixed_ht, &   not needed now

        V_East_FixedHeight, V_South_FixedHeight, V_Upward_FixedHeight, tts_fixed_ht, &   ! inputs

        !! telec_fixed_ht, &  not needed now

        inGIP, isGIP, &  ! was : IN, IS, &
        ! iwrite_plasma_interface, & not used lrm20121115
        TN_plasma_input_3d, &
        O_plasma_input_3d, &
        O2_plasma_input_3d, &
        N2_plasma_input_3d, &

        !! NO_plasma_input_3d, N4S_plasma_input_3d, N2D_plasma_input_3d, &  not needed now

        GLAt_plasma_3d, &
        GLOnd_plasma_3d, &
        PZ_plasma_3d, &
        V_east_plasma, V_south_plasma, V_upward_plasma, &  ! wind

        !! te_dum_plasma_input_3d, & not needed now

        ilon1_3d_fixed_ht, ilon2_3d_fixed_ht, &  ! output
        ilat1_3d_fixed_ht, ilat2_3d_fixed_ht, &  ! output
        ispecial_3d_fixed_ht, &
        ihl_3d_fixed_ht, ihu_3d_fixed_ht, &
        isFirstCallFixedHeight, &
        GIP_switches)

        CALL cpu_time(endTime)
        timeFixedGridtoFluxGrid = (endTime - startTime) + timeFixedGridtoFluxGrid 
        WRITE(*,*) "cpu_time for INTERFACE__FIXED_GEO_to_IONOSPHERE  : ",(endTime-startTime)

        !---------------------------------------
        ! print out results to check vs GT-IPE
        !---------------------------------------
        If (debugFixedGeo) then

            print *,'driver_ipe_gt.3d : AFTER CALLING INTERFACE__FIXED_GEO_to_MID_LAT_IONOSPHERE.............'

            CALL checkFixedGeo(debugdir, debugFixedGeoFileName, unitFixedGeo, &
                         npts, nmp, &
                         TN_plasma_input_3d, O_plasma_input_3d, O2_plasma_input_3d, &
                         N2_plasma_input_3d, GLAt_plasma_3d, GLOnd_plasma_3d, &
                         pz_plasma_3d, V_east_plasma, V_south_plasma, &
                         V_upward_plasma, ilon1_3d_fixed_ht, ilon2_3d_fixed_ht, &
                         ilat1_3d_fixed_ht, ilat2_3d_fixed_ht, ispecial_3d_fixed_ht, &
                         ihl_3d_fixed_ht, ihu_3d_fixed_ht, isFirstCallFixedHeight)

            print *,'STOPPING AFTER CHECKING FIXED GEO ******************************************'
            STOP

       endif ! debugFixedGeo


       !---------------------------------------
       ! Write out interpolated inputs to IPE
       !---------------------------------------
       if (debugThermoInterp) then 

        write(unitCheckThermoInterpFluxT,*) V_South_plasma
        write(unitCheckThermoInterpFluxT+1,*) V_East_plasma
        write(unitCheckThermoInterpFluxT+2,*) V_Upward_plasma
        write(unitCheckThermoInterpFluxT+3,*) ! should be rmt , but the composition is not being used ***
        write(unitCheckThermoInterpFluxT+4,*) TN_plasma_input_3d
        write(unitCheckThermoInterpFluxT+5,*) O_plasma_input_3d
        write(unitCheckThermoInterpFluxT+6,*) O2_plasma_input_3d
        write(unitCheckThermoInterpFluxT+7,*) N2_plasma_input_3d
        write(unitCheckThermoInterpFluxT+8,*) !qion3d_


       endif



else


    ! Update the neutral 3D structure: use MSIS/HWM to get the values in the flux tube grid
    CALL neutral ( utime )

end if  ! if ( sw_neutral == 'GT' )



   ! update plasma
   !t        CALL plasma ( utime )

   ! update electrodynamics
   !t        CALL eldyn ( utime )

   ! output to a file
   !t        CALL output ( utime )

   bigLoop = bigLoop + 1 !  increase # of times in big loop by 1

END DO  time_loop !: DO utime = start_time, stop_time, time_step

!-----------------------------------------------
! Average the times spent in the subroutines
!-----------------------------------------------
timeFluxGridtoFixedGrid = timeFluxGridtoFixedGrid/(bigLoop - 1)
timeFixedGridtoPressureGrid = timeFixedGridtoPressureGrid/(littleLoop - 1)
timePressureGridtoFixedGrid = timePressureGridtoFixedGrid/(bigLoop - 1)
timeFixedGridtoFluxGrid = timeFixedGridtoFluxGrid/(bigLoop - 1)

print *,'INTERPOLATION TIMES : -------------------------------------------------'
print *,'timeFluxGridtoFixedGrid = ', timeFluxGridtoFixedGrid
print *,'timeFixedGridtoPressureGrid = ', timeFixedGridtoPressureGrid
print *,'timePressureGridtoFixedGrid = ', timePressureGridtoFixedGrid
print *,'timeFixedGridtoFluxGrid = ', timeFixedGridtoFluxGrid
print *,' '
print *,'bigLoop, littleLoop = ', bigLoop, littleLoop


!-------------------------
!Close IPE startup files
!-------------------------
do ii = 1, numIonoStart
   
   CLOSE (UNIT=startUpFiles(ii)%unit)

enddo ! ii

!-------------------------------------------------------
! if we're in thermospheric debug mode, close the file
!-------------------------------------------------------
if (debugThermo) then
   CLOSE (unitCheckThermo)
endif

if (debugThermoInterp) CLOSE (unitCheckThermoInterp)

if (debugFixedGeo) CLOSE (unitFixedGeo)


! close all open files
!t      CALL close_files ()

print *,'driver_ipe_gt.3d : END OF driver_ipe_gt.3d '

END PROGRAM  test_GT
