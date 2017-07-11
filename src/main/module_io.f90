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
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: lun_min1=8000
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: lun_max1=lun_min1+ISPEC+3+ISPEV-1+1
      INTEGER (KIND=int_prec), DIMENSION(lun_min1:lun_max1),PUBLIC :: LUN_PLASMA1 !WRITE
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: lun_min2=9000
      INTEGER (KIND=int_prec), PARAMETER, PUBLIC :: lun_max2=lun_min2+ISPEC+3-1
      INTEGER (KIND=int_prec), DIMENSION(lun_min2:lun_max2),PUBLIC :: LUN_PLASMA2 !READ
      INTEGER (KIND=int_prec), PUBLIC :: record_number_plasma
!nm20120311
      INTEGER (KIND=int_prec), PUBLIC :: luntmp1,luntmp2,luntmp3
!nm20141001:
      INTEGER (KIND=int_prec), PUBLIC :: lun_wind0,lun_wind1,lun_wind2,lun_wind3,lun_wind4,lun_wind5

END MODULE module_IO
