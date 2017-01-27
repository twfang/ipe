!20111107: copied originally from GIP: ML__IONNEUT_PLAS
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
      MODULE module_IONNEUT_PLAS
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: IONNEUT_PLAS

      CONTAINS
!--------------------
      SUBROUTINE IONNEUT_PLAS(NPTS,P1,P2,P3,PI1,PI2,PI3,T,VIN,AMIn, &
     &  IN,IS,iout)
      USE module_precision
      IMPLICIT NONE
      INTEGER (KIND=int_prec), intent(in)  ::  NPTS
      REAL (KIND=real_prec) :: a , AMIn , amu , b , factor , P1 , P2 , P3 , PI1 , PI2 , &
           sum , summol , T , v1 , v2 , VIN , PI3
      INTEGER (KIND=int_prec) :: n , NMAx , iout, in, is
      
      DIMENSION P1(npts) , P2(npts) , P3(npts) , T(npts) , &
           VIN(npts) , AMIn(npts) , &
           a(3) , b(3) , PI1(npts) , PI2(npts) , PI3(npts)
      REAL (KIND=real_prec) :: mi1 , mi2 , mi3
      DATA mi1 , mi2 , mi3/16. , 30. , 32./
      DATA a/3.42E-11 , 6.66E-10 , 6.82E-10/
      DATA b/2.44E-10 , 4.28E-10 , 4.34E-10/
      amu = 1.66E-27
!c  **
!c  **
      factor=1.0
!c  **
!c  **
      DO 100 n = IN , IS
         summol = PI2(n) + PI3(n)
         sum = PI1(n) + PI2(n) + PI3(n)
         v2 = b(1)*P1(n) + b(2)*P2(n) + b(3)*P3(n)
         v1 = a(3)*P3(n) + a(2)*P2(n) + a(1)*P1(n)*factor*SQRT(T(n)) &
              *(1.08-0.139*log10(T(n))+4.51E-03*log10(T(n))**2)
         if(summol < 1.d-90) summol=0.0
         if(v1 < 1.d-90) v1=0.0
         if(v2 < 1.d-90) v2=0.0
         ! if(pi1(n).lt.1.d-90) pi1(n)=0.0
         ! if(iout.eq.1) write(6,*) 'here 5',n
         VIN(n) = (v1*PI1(n)+v2*summol)*1.E-06/sum
         AMIn(n) = (PI1(n)*mi1+PI2(n)*mi2+PI3(n)*mi3)*amu/sum
100   ENDDO
      RETURN

      END SUBROUTINE IONNEUT_PLAS
      END MODULE module_IONNEUT_PLAS
