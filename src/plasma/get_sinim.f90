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
      SUBROUTINE Get_sinI ( sw_sinI, sinI, GL_radi, D3 )
        USE module_precision
        IMPLICIT NONE
        INTEGER (KIND=int_prec),INTENT(IN)  :: sw_sinI
        REAL (KIND=real_prec),INTENT(IN)   :: GL_radi
        REAL (KIND=real_prec),INTENT(IN)   :: D3(3)
        REAL (KIND=real_prec),INTENT(OUT)  :: sinI
!---local variables---
        REAL (KIND=real_prec) :: CHX
        REAL (KIND=real_prec) :: SQTH
        REAL (KIND=real_prec) :: dotprod

!dbg20110831
!d print *,'sub-Get_sinI'
!d print *,'(2)',sw_sinI
!d print *,'(3)',sinI
!d print *,'(4)',GL_radi &
!d &,D3(1:3)

! (1) based on Phil's calculation
IF ( sw_sinI == 0 ) THEN
        CHX  = COS ( GL_radi )  !GL magnetic colatitude [rad]
        SQTH = SQRT( 3.0*(CHX*CHX)+1.0 )
        sinI  = 2.0*CHX/SQTH       !FD(2)

! (2) using APEX parameters
ELSE IF ( sw_sinI == 1 ) THEN 

!---
!!!DOT_PRODUCT( D3(1:3,i,mp), Vn_ms1(1:3,i) ) / SQRT(  DOT_PRODUCT( D3(1:3,i,mp), D3(1:3,i,mp) )  )
!---

        dotprod = D3(1)*D3(1) + D3(2)*D3(2) + D3(3)*D3(3)  !DOT_PRODUCT( D3,D3 )
        IF ( dotprod > 0.0 ) THEN
        sinI = D3(3) / &
     & SQRT(  dotprod   )
        ELSE
          sinI = 0.0
        END IF

END IF  !( sw_sinI == 0 ) THEN

      END SUBROUTINE Get_sinI
