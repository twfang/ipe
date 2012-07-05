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
      MODULE module_eldyn
      USE module_precision
      USE module_IPE_dimension,ONLY: NMP,NLP
!----------------------
!c idea
!      subroutine idea_geteb(im,ix,dayno,utsec,f107,kp,maglat,maglon,
!     &essa,ee1,ee2)
      USE efield
!c     use date_def
!c     use physcons, pi => con_pi
      IMPLICIT NONE
!      integer, intent(in) :: im  ! number of data points in efield 
!      integer, intent(in) :: ix  ! max data points in efield
!      integer :: dayno=254  ! calender day
!      real :: utsec=0.0  ! second
!      real :: f107=70.  ! 
!      real, intent(in) :: maglat(im)  ! magnetic latitude (rad)
!      real, intent(in) :: maglon(im)  ! magnetic longitude (rad)
!      real, intent(in) :: essa(im)  ! degree
!      real, intent(out)   :: ee1(im)    ! electric field x direction mV/m
!      real, intent(out)   :: ee2(im)    ! electric field y direction mV/m
!c     character*(*), intent(in) ::   dir    ! directory located coef files
!c local
!      integer i,k,iref,jref
!       real :: utsec_last
!      real utsec_last,dx,dy,aa,bb,maglond,maglatd,
!     &ed11(0:nmlon,0:nmlat),ed22(0:nmlon,0:nmlat),ylatm1(0:nmlat),
!     &ylonm1(0:nmlon),pi
!      real, parameter :: pi=3.141592653
!      logical first
!c
!      data first/.true./
!      data utsec_last/-1./
!      save first,utsec_last,ed11,ed22,ylatm1,ylonm1

      REAL(KIND=real_prec)   ,DIMENSION(0:nmlat  ),PUBLIC :: theta90_rad
      INTEGER(KIND=int_prec) ,DIMENSION(2  ,NLP  ),PUBLIC :: j0,j1 !1:NH; 2:SH
      REAL(KIND=real_prec)   ,DIMENSION(NMP,NLP*2),PUBLIC :: Ed1_90
      REAL(KIND=real_prec)   ,DIMENSION(NMP,NLP*2),PUBLIC :: Ed2_90
      REAL(KIND=real_prec)   ,DIMENSION(    NLP*2),PUBLIC :: coslam_m
      INTEGER (KIND=int_prec),DIMENSION(    NLP  ),PUBLIC :: lpconj !4 NH

!
      PRIVATE
      PUBLIC :: init_eldyn,eldyn
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
      SUBROUTINE eldyn ( utime )
      USE module_precision
      USE module_input_parameters,ONLY:NYEAR,NDAY,start_time
     &, ip_freq_output, sw_debug
     &, kp_eld
     &, F107D_ipe => F107D  !,AP
      USE module_physical_constants,ONLY:rtd
      IMPLICIT NONE
      INTEGER (KIND=int_prec),INTENT(IN)   :: utime !universal time [sec]
!---local
      real :: kp ! 
!c initiate
!c calculate efield only if diff time step
!      if(utsec.ne.utsec_last) then
!        utsec_last=utsec
      iday = NDAY !254 ! dayno                   ! day of year
!c       imo=idate(2)
!c       iday_m=idate(3) 
!!!need to create a routine to calculate month/day from NDAY!!!
      imo=3                     !month
      iday_m=15                 !day of month 
      iyear = NYEAR 
!!! F107D is global both in module efield & ipe input
      f107d = F107D_ipe         !f107
      ut = REAL(utime,real_prec)/3600.0
      kp = kp_eld  !=1.                   !???
      bz = .433726 - kp*(.0849999*kp + .0810363)  
     &        + f107d*(.00793738 - .00219316*kp)
      by=0.

      if ( utime==start_time ) then
        print *,'iday',iday, 'imo',imo,' iday_m',iday_m,' iyear',iyear
        print *,' utime=',utime,' kp',kp
        print *,'By=',by,' Bz=',bz,' F107d=',f107d
      end if

      call get_efield
!dbg        print*,'www'
!        ed11=ed1
!        ed22=ed2
      IF ( utime==start_time ) THEN
          write(unit=2003,FMT='(20f10.4)')ylatm
          write(unit=2004,FMT='(20f10.4)')ylonm
        END IF
!      endif !if(utsec.ne.utsec_last) then


! get ED1/2(nmp=80 X nlp=170) at 90km from potent(181x91)at 130km
        IF ( utime==start_time ) j0(:,:)=-999
        CALL GET_EFIELD90km ( utime )
        if ( sw_debug )  print *,'GET_EFIELD90km finished'
        IF ( utime==start_time ) write(unit=2007,FMT='(20f10.4)')
     &    (90.-theta90_rad*rtd)    

        IF ( MOD( (utime-start_time),ip_freq_output)==0 ) THEN
          write(unit=2000,FMT='(20E12.4)')potent !V
          write(unit=2001,FMT='(20E12.4)')ed1 *1.0E+03 !V/m-->mV/m
          write(unit=2002,FMT='(20E12.4)')ed2 *1.0E+03 !V/m-->mV/m
          write(unit=2008,FMT='(20E12.4)')ed1_90 *1.0E+03 !V/m-->mV/m
          write(unit=2009,FMT='(20E12.4)')ed2_90 *1.0E+03 !V/m-->mV/m
          write(unit=2010,FMT='(I12)')utime !sec
        END IF
!c
      return
!      end
!      end program ts_efield
!
      END SUBROUTINE eldyn
      END MODULE module_eldyn
