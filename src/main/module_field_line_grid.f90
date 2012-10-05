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
!SMS$DISTRIBUTE(dh,2,1) BEGIN
      INTEGER(KIND=int_prec),ALLOCATABLE :: JMIN_IN_all(:,:),JMAX_IS_all(:,:)!.. first and last indices on field line grid
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,1) BEGIN
      INTEGER(KIND=int_prec),ALLOCATABLE,PUBLIC :: JMIN_IN (:),JMAX_IS (:)   !.. first and last indices on field line grid
      INTEGER(KIND=int_prec),ALLOCATABLE,PUBLIC :: JMIN_ING(:),JMAX_ISG(:)   !.. first and last indices on field line grid
!SMS$DISTRIBUTE END
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
      INTEGER (KIND=int_prec) :: ISL=1,IBM=2,IGR=3,IQ=4,IGCOLAT=5,IGLON=6
!SMS$DISTRIBUTE(dh,2,3) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: plasma_grid_3d(:,:,:,:)      !(MaxFluxTube,NLP,NMP,6)
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: Be3(:,:,:) ! .. Eq(4.13) Richmond 1995 at Hr=90km in the NH(1)/SH(2) foot point [T]
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,dimension(:,:,:) :: ON_m3,HN_m3,N2N_m3,O2N_m3,HE_m3,N4S_m3,TN_k,TINF_k
!SMS$DISTRIBUTE END
!V_ExB m/s at the apex height
!SMS$DISTRIBUTE(dh,1,2) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,DIMENSION(:,:) :: VEXBup !DIMENSION(NMP,NLP)
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,3,4) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: plasma_3d(:,:,:,:) !ISTOT,MaxFluxTube,NLP,NMP
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: Un_ms1   (:,:,:,:) !Ue1 Eq(5.6) in magnetic frame !1st dim: corresponds to apexD1-3
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: apexD  (:,:,:,:,:) !(3,MaxFluxTube,NLP,NMP,3).. Eq(3.8-10 ) Richmond 1995
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: apexE  (:,:,:,:,:) !(2,MaxFluxTube,NLP,NMP,3).. Eq(3.11-12) Richmond 1995
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,2) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: plasma_grid_Z (:,:)  !.. altitude [meter] (MaxFluxTube,NLP)
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: plasma_grid_GL(:,:)  !.. magnetic co-latitude Eq(6.1) [rad]
      REAL(KIND=real_prec),ALLOCATABLE               :: r_meter2D     (:,:)  !.. distance from the center of the Earth[meter]
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,3) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: hrate_cgs_save(:,:,:)!.. each component of the Neutral heating rate (eV/cm^3/s) DIM(7,MaxFluxTube,NMP)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,,1) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: mlon_rad(:)          !mag longitude in [rad]
!SMS$DISTRIBUTE END

      REAL(KIND=real_prec),PARAMETER  ,PUBLIC        :: dlonm90km = 4.5 ![deg]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  GCOLAT(:,:)    !.. geographic co-latitude [rad]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  GLON(:,:)      !.. geographic longitude [rad]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  Qvalue(:,:)
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  GL_rad(:,:)    !.. magnetic co-latitude Eq(6.1) [rad]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  SL_meter(:,:)  !.. distance of point from northern hemisphere foot point [meter]
!      REAL(KIND=real_prec), ALLOCATABLE, PUBLIC ::  BM_T(:,:)      !.. magnetic field strength [T]

!------------
!...calculated parameters
!      REAL(KIND=real_prec), ALLOCATABLE,     PUBLIC ::  Z_meter(:,:) !.. altitude [meter]
!SMS$DISTRIBUTE(dh,1) BEGIN
      REAL(KIND=real_prec), ALLOCATABLE,     PUBLIC ::  Pvalue(:)  !.. p coordinate (L-shell)
!SMS$DISTRIBUTE END
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
!
!nm20121003:subroutines read_plasma_grid_global init_plasma_grid are separated into module_sub_field_line_grid.f90.


!---------------------------
      END MODULE module_FIELD_LINE_GRID_MKS
