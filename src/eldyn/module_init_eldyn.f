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
!Aug2011: the original code was provided from Fei Wu from WAM version
!nm20110906: modified to implement to IPE
!copied from /.../testall.f
!ylonm(1:nmlon=180)
!ylatm(1:nmlat=90)
!      program ts_efield
      MODULE module_init_eldyn
      USE module_precision
      USE module_IPE_dimension,ONLY: NMP,NLP
!----------------------
!c idea
!      subroutine idea_geteb(im,ix,dayno,utsec,f107,kp,maglat,maglon,
!     &essa,ee1,ee2)
      USE efield !,ONLY:iday,imo,iday_m,iyear,ut,kp,by,bz,f107d
      USE module_efield_init,ONLY:efield_init
!c     use date_def
!c     use physcons, pi => con_pi
      IMPLICIT NONE

!nm20121003:module parameters are separated into module_eldyn.f

      PRIVATE
      PUBLIC :: init_eldyn
      CONTAINS
      SUBROUTINE init_eldyn ( )
      IMPLICIT NONE
!20120304:      CHARACTER(len=*),PARAMETER :: path='~/sandbox/efield/'
!
      print *,'begin init_eldyn'
      CALL efield_init( 
     &'coeff_lflux.dat',
     &'coeff_hflux.dat',
     &'wei96.cofcnts')
!
      print *,'END sub-init_eld'
      END SUBROUTINE init_eldyn
!---
      END MODULE module_init_eldyn
