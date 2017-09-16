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
      MODULE module_sub_eldyn
      use module_precision
      USE efield !,ONLY:iday,imo,iday_m,iyear,ut,kp,by,bz,f107d
      use module_get_efield,ONLY:get_efield
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: eldyn
      CONTAINS

      SUBROUTINE eldyn ( utime )
      use module_precision
      use module_cal_monthday
      use module_input_parameters,ONLY:NYEAR,NDAY,start_time,mype       &
     &,ip_freq_output,sw_debug,kp_eld,sw_eldyn,F107D_ipe => F107D       !,AP
      use module_physical_constants,ONLY:rtd
      use module_eldyn,ONLY:theta90_rad,j0,Ed1_90,Ed2_90
      use module_sunloc,only: sunloc
      use module_highlat,only: highlat
      use module_sub_dynamo,only: dynamo
      use module_magfield,ONLY:sunlons
      use module_update_fli,ONLY:update_fli
      use dynamo_module,only:zigm11

      IMPLICIT NONE
      INTEGER (KIND=int_prec),INTENT(IN)   :: utime !universal time [sec]
!---local
      real :: kp ! 
      integer (KIND=int_prec) :: iyr
      real (KIND=real_prec)   :: utsecs
      character :: fname*10,labl*56,units*12

      print *,'begin sub_eldyn: sw_eldyn=', sw_eldyn
!1: self-consistent electrodynamic solver
!t      IF ( sw_eldyn==0 ) THEN 

!t      n_time=n_time+1
      print *, 'self-consistent eldyn started' !t ,n_time

      iyr = 1997
      utsecs=REAL(utime, real_prec)

      
!t output time in netcdf
!t      fname = 'time'
!t      labl = 'time'
!t      units = 'sec'
!t      call ncplot1D(noid,fname,labl,start1t,count1t,idtime              &
!t     &     ,utsecs,4,units,3)



      print *,'sub-eldyn: sunloc: utsecs=',iyr,NDAY,utsecs
!      call sunloc(iyr,NDAY,utsecs)
      call sunloc(iyr,97,utsecs)

!dbg20150615: temporary commented out
! output sunlons(1)
!      fname = 'sunlons'
!      labl = 'sunlons'
!      units = 'radian'
!      start1_out(1)=n_time
!      count1(1)    =n_time
!      dim1         =dim3(3)
!      call ncplot1D(noid,fname,labl,start1_out,count1,dim1              &
!     &     ,sunlons,7,units,6)

! output sunlons
      print *,'(21) output dyn sunlons at utime=',utime
!      write(unit=4021,FMT='(I12)')utime
      write(unit=4021,FMT='(20E12.4)')sunlons(1)
      if (sw_debug) print *,'sunlons(1)',sunlons(1)

      if (sw_debug) print *,'sub-eldyn: update_fli'
      call update_fli ( utime )

      if (sw_debug) print *,'sub-eldyn: highlat'
      call highlat

      if (sw_debug) print *,'sub-dynamo: dynamo'
      call dynamo
      print *,'self-consistent dynamo finished'

!2: WACCM empirical electric field model
!t      ELSE IF ( sw_eldyn==1 ) THEN 
!c initiate
!c calculate efield only if diff time step
!      if(utsec.ne.utsec_last) then
!        utsec_last=utsec
      iday = NDAY !254 ! dayno                   ! day of year
!c       imo=idate(2)
!c       iday_m=idate(3) 
!!!need to create a routine to calculate month/day from NDAY!!!
!      imo=3                     !month
      iyear = NYEAR 
!nm20121127: calculate month/day from iyear and iday
      call cal_monthday ( iyear,iday, imo,iday_m )
!nm20130402: temporarily hard-code the iday to get b4bconfirmed.
      iday_m=15                 !day of month 

!!! F107D is global both in module efield & ipe input
      f107d = F107D_ipe         !f107
      ut = REAL(utime,real_prec)/3600.0 !sec-->hr
      kp = kp_eld  !=1.                   !???
      bz = .433726 - kp*(.0849999*kp + .0810363)                        &
     &        + f107d*(.00793738 - .00219316*kp)
      by=0.

      if ( utime==start_time ) then
        print *,'iday',iday, 'imo',imo,' iday_m',iday_m,' iyear',iyear
        print *,' kp',kp
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
      IF ( utime==start_time ) then
        j0=-999
      endif
      CALL GET_EFIELD90km ( utime )
      if ( sw_debug )  print *,'GET_EFIELD90km finished'

!t      END IF !( sw_eldyn==0 ) THEN 


      IF ( utime==start_time ) THEN 
        write(unit=2007,FMT='(20f10.4)') (90.-theta90_rad*rtd)    
      ENDIF

      IF ( MOD( (utime-start_time),ip_freq_output)==0 ) THEN
!SMS$SERIAL(<ed1_90,ed2_90,IN>:default=ignore) BEGIN
        write(unit=2000,FMT='(20E12.4)')potent !V
        write(unit=2001,FMT='(20E12.4)')ed1 *1.0E+03 !V/m-->mV/m
        write(unit=2002,FMT='(20E12.4)')ed2 *1.0E+03 !V/m-->mV/m
        write(unit=2008,FMT='(20E12.4)')ed1_90 *1.0E+03 !V/m-->mV/m
        write(unit=2009,FMT='(20E12.4)')ed2_90 *1.0E+03 !V/m-->mV/m
        write(unit=2010,FMT='(I12)')utime !sec
 !SMS$SERIAL END
       END IF

      return

      END SUBROUTINE eldyn
      END MODULE module_sub_eldyn
