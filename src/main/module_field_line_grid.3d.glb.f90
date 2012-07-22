!dbg2010923: to reduce memory1: apexE etc
!dbg20110927: to reduce memory2:
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
      MODULE module_FIELD_LINE_GRID_MKS
      USE module_precision
      USE module_IPE_dimension,ONLY: NMP,NLP,ISTOT
      IMPLICIT NONE

! --- PRIVATE ---
!
! --- PUBLIC ---
      INTEGER (KIND=int_prec),PUBLIC :: mp_save,lp_save,MaxFluxTube
      REAL(KIND=real_prec),PARAMETER,PUBLIC :: ht90  = 90.0E+03  !reference height in meter
!... read in parameters
      INTEGER(KIND=int_prec), PUBLIC,DIMENSION(NMP,NLP) :: JMIN_IN_all,JMAX_IS_all  !.. first and last indices on field line grid
      INTEGER(KIND=int_prec), ALLOCATABLE,TARGET,PUBLIC :: JMIN_IN(:),JMAX_IS(:)    !.. first and last indices on field line grid
      TYPE :: plasma_grid
!dbg20110927         REAL(KIND=real_prec) :: Z  !.. altitude [meter]
         REAL(KIND=real_prec) :: SL !.. distance of point from northern hemisphere foot point [meter]
         REAL(KIND=real_prec) :: BM !.. magnetic field strength [T]
         REAL(KIND=real_prec) :: GR !.. Gravity [m2 s-1]
!dbg20110927         REAL(KIND=real_prec) :: GL !.. magnetic co-latitude Eq(6.1) [rad]
         REAL(KIND=real_prec) :: Q  
         REAL(KIND=real_prec) :: GCOLAT !.. geographic co-latitude [rad]
         REAL(KIND=real_prec) :: GLON   !.. geographic longitude [rad]
      END TYPE plasma_grid
!     TYPE(plasma_grid), ALLOCATABLE, TARGET,PUBLIC :: plasma_grid_3d(:,:)      
      INTEGER (KIND=int_prec) :: ISL=1,IBM=2,IGR=3,IQ=4,IGCOLAT=5,IGLON=6
!SMS$DISTRIBUTE(dh,1,2) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,TARGET,PUBLIC :: plasma_grid_3d(:,:,:,:)      
!SMS$DISTRIBUTE END
!V_ExB m/s at the apex height
      REAL(KIND=real_prec), DIMENSION(:,:),ALLOCATABLE, PUBLIC :: VEXBup !DIMENSION(NMP,NLP)

!SMS$DISTRIBUTE(dh,1) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,TARGET,PUBLIC :: plasma_3d(:,:,:,:)   !ISTOT,MaxFluxTube,NLP,NMP
      REAL(KIND=real_prec),ALLOCATABLE,TARGET,PUBLIC :: plasma_grid_Z (:,:)  !.. altitude [meter] (MaxFluxTube,NLP)
      REAL(KIND=real_prec),ALLOCATABLE,TARGET,PUBLIC :: plasma_grid_GL(:,:)  !.. magnetic co-latitude Eq(6.1) [rad]
      REAL(KIND=real_prec),ALLOCATABLE               :: r_meter2D     (:,:)  !.. distance from the center of the Earth[meter]
      REAL(KIND=real_prec),ALLOCATABLE,TARGET,PUBLIC :: mlon_rad(:)          !mag longitude in [rad]
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: hrate_cgs_save(:,:,:)!.. each component of the Neutral heating rate (eV/cm^3/s) DIM(7,NPTS2D,NMP)
!SMS$DISTRIBUTE END
      REAL(KIND=real_prec),allocatable,dimension(:,:,:),PUBLIC:: ON_m3,HN_m3,N2N_m3,O2N_m3,HE_m3,N4S_m3,TN_k,TINF_k
      REAL(KIND=real_prec),allocatable,dimension(:,:,:,:),PUBLIC  :: Un_ms1  !Ue1 Eq(5.6) in magnetic frame !1st dim: corresponds to apexD1-3

      REAL(KIND=real_prec),PARAMETER, PUBLIC :: dlonm90km = 4.5 ![deg]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  GCOLAT(:,:)    !.. geographic co-latitude [rad]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  GLON(:,:)      !.. geographic longitude [rad]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  Qvalue(:,:)
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  GL_rad(:,:)    !.. magnetic co-latitude Eq(6.1) [rad]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  SL_meter(:,:)  !.. distance of point from northern hemisphere foot point [meter]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  BM_T(:,:)      !.. magnetic field strength [T]

!SMS$DISTRIBUTE(dh,2,3) BEGIN
      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  Be3(:,:,:)     ! .. Eq(4.13) Richmond 1995 at Hr=90km in the NH(1)/SH(2) foot point [T]
!SMS$DISTRIBUTE END
!------------
!...calculated parameters
!      REAL(KIND=real_prec), ALLOCATABLE,     PUBLIC ::  Z_meter(:,:) !.. altitude [meter]
      REAL(KIND=real_prec), ALLOCATABLE,     PUBLIC ::  Pvalue(:)  !.. p coordinate (L-shell)
!      REAL(KIND=real_prec), ALLOCATABLE,     PUBLIC ::  GR_mks(:,:)  !.. Gravity [m2 s-1]


!-------------
!nm20110822:no longer used
!      REAL(KIND=real_prec),               ALLOCATABLE, PUBLIC ::  SZA_rad(:) !solar zenith angle [radians]
!
! components (east, north, up) of base vectors
      TYPE :: geographic_coords
         REAL(KIND=real_prec) :: east
         REAL(KIND=real_prec) :: north
         REAL(KIND=real_prec) :: up
      END TYPE geographic_coords

      INTEGER (KIND=int_prec) :: east=1,north=2,up=3
!SMS$DISTRIBUTE(dh,2,3) BEGIN
      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  apexD(:,:,:,:,:) !(3,MaxFluxTube,NLP,NMP,3).. Eq(3.8-10 ) Richmond 1995
      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  apexE(:,:,:,:,:) !(2,MaxFluxTube,NLP,NMP,3).. Eq(3.11-12) Richmond 1995
!SMS$DISTRIBUTE END
!
      PRIVATE :: read_plasma_grid_global !dbg ,test_grid_output
      PUBLIC :: init_plasma_grid


      CONTAINS
!---------------------------
! initialise plasma grids
        SUBROUTINE init_plasma_grid ( )
        USE module_physical_constants,ONLY: earth_radius, pi, G0,zero
        USE module_input_parameters,ONLY: sw_debug,lpstrt,lpstop,mpstrt,mpstop,sw_grid  
!dbg, fac_BM
        IMPLICIT NONE

        INTEGER (KIND=int_prec) :: i, mp,lp
        REAL (KIND=real_prec) :: sinI
        INTEGER (KIND=int_prec), parameter :: sw_sinI=0  !0:flip; 1:APEX
        INTEGER (KIND=int_prec) :: in,is
        INTEGER (KIND=int_prec) :: midpoint
!---------
! array initialization

!if the new GLOBAL 3D version
      CALL read_plasma_grid_global

! make sure to use the MKS units.
      print *,"Z_meter calculation completed"

      Pvalue(:) = zero
!SMS$PARALLEL(dh, mp, lp) BEGIN
      apex_longitude_loop: DO mp = 1,NMP
      IF ( sw_debug.AND.mpstrt<=mp.AND.mp<=mpstop ) & 
     &  print *,"sub-init_plasma_grid: mp=",mp


!.. p coordinate (L-shell) is a single value along a flux tube 
!NOTE: in FLIP, PCO is only used in setting up the rough plasmasphere H+ initial profiles (See PROFIN). It does not have to be accurate.

!dbg20120112:      Pvalue(:) = zero
      apex_latitude_height_loop:   DO lp = 1,NLP

        IN = JMIN_IN(lp)
        IS = JMAX_IS(lp)

        IF (mp==1)  CALL Get_Pvalue_Dipole ( r_meter2D(IN,lp), plasma_grid_GL(IN,lp), Pvalue(lp) )

!debug write
IF ( sw_debug.AND. & !) THEN
& mpstrt<=mp.AND.mp<=mpstop .AND. lpstrt<=lp.AND.lp<=lpstop) THEN

!dbg20120305
midpoint = IN + ( IS - IN )/2
!print *,'midpoint',midpoint,plasma_grid_Z(midpoint)

print "('lp=',i6,'  IN=',i6,'  IS=',i6,'  NPTS=',i6)", lp,IN,IS,(IS-IN+1)
print "('r [m]      =',2E12.4)", r_meter2D(in,lp),r_meter2D(is,lp)
print "('G-LAT [deg]=',2f10.4)",(90.-plasma_grid_3d(in,lp,mp,IGCOLAT)*180./pi),(90.-plasma_grid_3d(is,lp,mp,IGCOLAT)*180./pi)
print "('M-LAT [deg]=',2f10.4)",(90.-plasma_grid_GL(in,lp)*180./pi),(90.-plasma_grid_GL(is,lp)*180./pi)
print "('GLON  [deg]=',2f10.4)",(plasma_grid_3d(in,lp,mp,IGLON)*180./pi),(plasma_grid_3d(is,lp,mp,IGLON)*180./pi)
print "('Qvalue     =',2E12.4)", plasma_grid_3d(in,lp,mp,IQ), plasma_grid_3d(is,lp,mp,IQ)
print "('BM [Tesla]    =',2E12.4)", plasma_grid_3d(in,lp,mp,IBM), plasma_grid_3d(is,lp,mp,IBM)
!print "('D1         =',2E12.4)", apexD(1,in,mp)%east, apexD(1,is,mp)%east
!print "('D2         =',2E12.4)", apexD(2,in,mp)%east, apexD(2,is,mp)%east
print "('D3         =',6E12.4)", apexD(3,in,lp,mp,east),apexD(3,in,lp,mp,north),apexD(3,in,lp,mp,up), apexD(3,is,lp,mp,east),apexD(3,is,lp,mp,north),apexD(3,is,lp,mp,up)
print "('E1         =',2E12.4)", apexE(1,in,lp,mp,east), apexE(1,is,lp,mp,east)
print "('E2         =',2E12.4)", apexE(2,in,lp,mp,east), apexE(2,is,lp,mp,east)
print "('Be3 [T] NH/SH  =',2E12.4)", Be3(1,mp,lp), Be3(2,mp,lp)

print "('SL [m]     =',4E13.5)", plasma_grid_3d(in:in+1,lp,mp,ISL), plasma_grid_3d(is-1:is,lp,mp,ISL)
print "('Z  [m]     =',4E13.5)",  plasma_grid_Z(in:in+1,lp),  plasma_grid_Z(is-1:is,lp)
print "('Pvalue     =',F10.4)", Pvalue(lp)
END IF !( sw_debug etc ) THEN


IF ( sw_grid==0 ) THEN  !APEX
! assuming Newtonian gravity: G0 is gravity at the sea level (z=0) 
!NOTE: positive in NORTHern hemisphere; negative in SOUTHern hemisphere
         flux_tube: DO i=IN,IS

!dbg20110831
!d print *,'calling sub-Get_sinI'&
!d &, i,mp,lp,sw_sinI, sinI&
!d &, plasma_grid_GL(i), apexD(3,i,mp)%east, apexD(3,i,mp)%north, apexD(3,i,mp)%up

           CALL Get_sinI ( sw_sinI, sinI, plasma_grid_GL(i,lp) &
     &, apexD(3,i,lp,mp,east), apexD(3,i,lp,mp,north), apexD(3,i,lp,mp,up) ) 
           plasma_grid_3d(i,lp,mp,IGR)  =  G0 * ( earth_radius * earth_radius ) / ( r_meter2D(i,lp) * r_meter2D(i,lp) ) * sinI * (-1.0)

!IF ( sw_debug )  print "(4E12.4)", Z_meter(i,mp),sinI, (G0 * ( earth_radius * earth_radius ) / ( r_meter2D(i,lp) * r_meter2D(i,lp) )),  GR_mks(i,mp)
        
         END DO flux_tube

IF ( sw_debug )  then
  if ( sw_sinI==0 ) then
    print *,'sinI: flip'
  else if ( sw_sinI==1 ) then
    print *, 'sinI: APEX'
  endif !  if ( sw_sinI==0 ) then
  print "('GRavity[m2 s-1]=',4E12.4)",plasma_grid_3d(in:in+2,lp,mp,IGR),plasma_grid_3d(is,lp,mp,IGR)
END IF !( sw_debug )  then

!nm20120304: introducing the flip grid
!reminder:
!(1) neutral wind should be calculated using sinI from flip_grid: "SINDIP"
ELSE IF ( sw_grid==1 ) THEN  !FLIP
  IF ( mpstrt<=mp.AND.mp<=mpstop .AND. lpstrt<=lp.AND.lp<=lpstop) THEN

   print *,'calling get_FLIP_grid',mp,lp
   CALL get_flip_grid (mp,lp)
  ENDIF
END IF !( sw_grid==0 ) THEN  !APEX

!debug20110314
if( sw_debug .and. mp==1 .and. lp>=lpstrt .and. lp<=lpstop ) then
print *,'lp=',lp,' in=',in,' is=',is,(is-in+1),(90.-plasma_grid_GL(in,lp)*180./pi)
endif !(mp==1) then

       END DO apex_latitude_height_loop   !: DO lp = 1,NLP
     END DO apex_longitude_loop         !: DO mp = 1,NMP 
!SMS$PARALLEL END

     mlon_rad(:) = zero
     DO mp = 1,NMP+1 
       mlon_rad(mp) = REAL( (mp-1),real_prec ) * dlonm90km *pi/180.00
     END DO
if ( sw_debug ) print *,'mlon_rad[deg]',mlon_rad*180.0/pi

!dbg20120313
!     DO mp = 1,NMP 
!       DO i = 1,NPTS2D 
!         plasma_grid_3d(i,lp,mp,IBM) = plasma_grid_3d(i,lp,mp,IBM) * fac_BM
!       END DO
!     END DO

        END SUBROUTINE init_plasma_grid
!---------------------------
!20110726: the new GLOBAL 3D version: NMP=80
! global grid with the low resolution version
! magnetic longitude used for the grid. from 0 to 355.5 with 4.5 degree interval

      SUBROUTINE read_plasma_grid_global
        USE module_IPE_dimension,ONLY: NPTS2D,NMP,NLP
        USE module_physical_constants,ONLY: earth_radius,pi,zero
        USE module_input_parameters,ONLY:read_input_parameters,sw_debug,lpstrt,sw_neutral_heating_flip
        USE module_IO,ONLY: filename,LUN_pgrid
        IMPLICIT NONE

!        integer (kind=int_prec), parameter :: NMP=80
!        integer (kind=int_prec), parameter :: NLP=170
!        integer (kind=int_prec), parameter :: NPTS2D=44438 
!        real(kind=8) gr_2d(npts2,nmp)
!        real(kind=8) gcol_2d(npts2,nmp)
!    real(kind=8) glon_2d(npts2,nmp)
!    real(kind=8) q_coordinate_2d(npts2,nmp)
!    real(kind=8) bcol_2d(npts2,nmp)
!    real(kind=8) Apex_D1_2d(3,npts2,nmp)
!    real(kind=8) Apex_D2_2d(3,npts2,nmp)
!    real(kind=8) Apex_D3_2d(3,npts2,nmp)
!    real(kind=8) Apex_E1_2d(3,npts2,nmp)
!    real(kind=8) Apex_E2_2d(3,npts2,nmp)
!    real(kind=8) Apex_grdlbm2_2d(3,npts2,nmp)
!    real(kind=8) integral_ds_2d(npts2,nmp)
!    real(kind=8) apex_BMAG_2d(npts2,nmp)
!    real(kind=8)  Apex_BE3_N(nmp,nlp), Apex_BE3_S(nmp,nlp)

!-------------
!... read in parameters
      INTEGER(KIND=int_prec) lp,mp,stat_alloc

      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  dum0    !.. distance from the center of the Earth[meter]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  dum1    !.. geographic co-latitude [rad]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  dum2    !.. geographic longitude [rad]
      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  dum3
!dbg20110927      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  GL_rad_all    !.. magnetic co-latitude Eq(6.1) [rad]
!dbg20110927      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  SL_meter_all  !.. distance of point from northern hemisphere foot point [meter]
!dbg20110927      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP) ::  BM_T_all      !.. magnetic field strength [T]
! components (east, north, up) of base vectors
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  dum4    !.. Eq(3.8) Richmond 1995
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  dum5    !.. Eq(3.9) Richmond 1995
      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  dum6    !.. Eq(3.10) Richmond 1995
!dbg20110927      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  E1_all    !.. Eq(3.11) Richmond 1995
!dbg20110927      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP) ::  E2_all    !.. Eq(3.12) Richmond 1995
!JFM  REAL(KIND=real_prec), DIMENSION(2,NMP,NLP) ::  Be3_all         ! .. Eq(4.13) Richmond 1995 at Hr=90km in the NH(1)/SH(2) foot point [T]
      REAL(KIND=real_prec), DIMENSION(NMP,NLP) ::  Be3_all1,Be3_all2 ! .. Eq(4.13) Richmond 1995 at Hr=90km in the NH(1)/SH(2) foot point [T]

!-------------local
        CHARACTER (LEN=11) :: FORM_dum
        CHARACTER (LEN=7)  :: STATUS_dum
        CHARACTER(LEN=*), PARAMETER :: filepath_pgrid= &
!     & '../../field_line_grid/20110419lowres_global/'  !20110419 low res global grid
     &  './'
        CHARACTER(LEN=*), PARAMETER :: filename_pgrid= &
!     &  'GIP_apex_coords_global_lowresCORRECTED'        !20110824: corrected for lp=1-6
     &  'ipe_grid'        !20110824: corrected for lp=1-6

!dbg      INTEGER(KIND=int_prec) :: sw_test_grid=0  !1: ON testgrid; 0: OFF 
!---

!SMS$SERIAL BEGIN
      filename =filepath_pgrid//filename_pgrid
      FORM_dum ='formatted' 
      STATUS_dum ='old'
      CALL open_file ( filename, LUN_pgrid, FORM_dum, STATUS_dum ) 
 
      print *,"open file completed"
      READ (UNIT=LUN_pgrid, FMT=*) JMIN_IN_all, JMAX_IS_all  !IN_2d_3d , IS_2d_3d
      print *,"reading JMIN_IN etc completed"
!SMS$SERIAL END

      MaxFluxTube = maxval(JMAX_IS_all(1,:)-JMIN_IN_all(1,:)+1)

      CALL allocate_arrays ( 0 )

      JMIN_IN = 1
      JMAX_IS = JMAX_IS_all(1,:) - JMIN_IN_all(1,:) + 1

! array initialization
      Be3            = zero
      plasma_grid_3d = zero
      plasma_grid_Z  = zero
      plasma_grid_GL = zero
      plasma_3d      = zero
      apexD          = zero
      apexE          = zero
      r_meter2D      = zero

!SMS$SERIAL(<dum0,dum1,dum2,dum3,r_meter2D,plasma_grid_3d,plasma_grid_GL,OUT>) BEGIN
      READ (UNIT=LUN_pgrid, FMT=*) dum0, dum1, dum2, dum3 !gr_2d, gcol_2d, glon_2d, q_coordinate_2d
do lp=1,NLP
  r_meter2D    (JMIN_IN(lp):JMAX_IS(lp),lp) = dum0(JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1)                !r_meter
  plasma_grid_Z(JMIN_IN(lp):JMAX_IS(lp),lp) = dum0(JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1) - earth_radius ![meter]
enddo
do mp=1,NMP
  do lp=1,NLP
    plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,mp,IGCOLAT) = dum1(JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),mp) !GCOLAT
    plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,mp,IGLON  ) = dum2(JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),mp) !GLON
    plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,mp,IQ     ) = dum3(JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),mp) !Q
  enddo
enddo
      print *,"reading r_meter etc completed"

      READ (UNIT=LUN_pgrid, FMT=*) dum0          !bcol_2d
do lp=1,NLP
  plasma_grid_GL(JMIN_IN(lp):JMAX_IS(lp),lp) = dum0(JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1) !GL
enddo
      print *,"reading GL_rad etc completed"
      READ (UNIT=LUN_pgrid, FMT=*) dum0, dum1 !integral_ds_2d, apex_BMAG_2d
do lp=1,NLP
  plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,ISL) = dum0(JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP) !SL
  plasma_grid_3d(JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,IBM) = dum1(JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP) !BM
enddo
      print *,"reading SL_meter etc completed"
!SMS$SERIAL END
!SMS$SERIAL(<dum4,dum5,dum6,apexD,OUT>) BEGIN
      READ (UNIT=LUN_pgrid, FMT=*) dum4, dum5, dum6      !Apex_D1_2d
!D2
!dbg20110923  apexD(1,1:NPTS2D,1:NMP)%east  =  dum4(1,1:NPTS2D,1:NMP) !D1
!dbg20110923  apexD(1,1:NPTS2D,1:NMP)%north =  dum4(2,1:NPTS2D,1:NMP)
!dbg20110923  apexD(1,1:NPTS2D,1:NMP)%up    =  dum4(3,1:NPTS2D,1:NMP)
!D2
!dbg20110923  apexD(2,1:NPTS2D,1:NMP)%east  =  dum5(1,1:NPTS2D,1:NMP) !D2
!dbg20110923  apexD(2,1:NPTS2D,1:NMP)%north =  dum5(2,1:NPTS2D,1:NMP)
!dbg20110923  apexD(2,1:NPTS2D,1:NMP)%up    =  dum5(3,1:NPTS2D,1:NMP)
!D3
do lp=1,NLP
  apexD(3,JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,east ) =  dum6(1,JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP) !D3
  apexD(3,JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,north) =  dum6(2,JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP)
  apexD(3,JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,up   ) =  dum6(3,JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP)
enddo
      print *,"reading D1-3 etc completed"
!SMS$SERIAL END
!SMS$SERIAL(<dum4,dum5,apexE,OUT>) BEGIN
      READ (UNIT=LUN_pgrid, FMT=*) dum4, dum5          !Apex_E1_2d
!E1
do lp=1,NLP
  apexE(1,JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,east ) =  dum4(1,JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP) !E1
  apexE(1,JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,north) =  dum4(2,JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP)
  apexE(1,JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,up   ) =  dum4(3,JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP)
!E2
  apexE(2,JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,east ) =  dum5(1,JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP) !E2
  apexE(2,JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,north) =  dum5(2,JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP)
  apexE(2,JMIN_IN(lp):JMAX_IS(lp),lp,1:NMP,up   ) =  dum5(3,JMIN_IN_all(1,lp):JMAX_IS_all(1,lp),1:NMP)
enddo
      print *,"reading E1/2 etc completed"
!SMS$SERIAL END
!SMS$SERIAL(<Be3,OUT>) BEGIN
!JFM  READ (UNIT=LUN_pgrid, FMT=*) Be3_all(1,1:NMP,1:NLP),Be3_all(2,1:NMP,1:NLP) !Apex_BE3_N
!JFM  READ (UNIT=LUN_pgrid, FMT=*) Be3_all(1,:,:),Be3_all(2,:,:) !Apex_BE3_N
      READ (UNIT=LUN_pgrid, FMT=*) Be3_all1,Be3_all2 !Apex_BE3_N
!JFM  Be3(1:2,1:NMP,1:NLP)=    Be3_all(1:2,1:NMP,1:NLP)
Be3(1,1:NMP,1:NLP)=    Be3_all1(1:NMP,1:NLP)
Be3(2,1:NMP,1:NLP)=    Be3_all2(1:NMP,1:NLP)
      print *,"reading Be3 etc completed"
      CLOSE(UNIT=LUN_pgrid)
      print *,"global grid reading finished, file closed..."
!SMS$SERIAL END


!dbg20110811:
!dbg IF ( sw_test_grid==1 ) THEN

!dbg  CALL test_grid_output ( JMIN_IN_all,JMAX_IS_all,r_meter_all,SL_meter_all,BM_T_all,GL_rad_all &
!dbg & ,GCOLAT_all,GLON_all,Qvalue_all,D1_all,D2_all,D3_all,E1_all,E2_all,Be3_all &
!dbg &,filepath_pgrid,filename_pgrid)
!dbg END IF 

      END SUBROUTINE  read_plasma_grid_global
!---
!20110927:test_grid_output deleted
!---------------------------
      END MODULE module_FIELD_LINE_GRID_MKS
