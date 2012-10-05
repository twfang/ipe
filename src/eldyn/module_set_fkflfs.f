!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_set_fkflfs
!--------------------------------------------------------------------- 
! description: calculates the electric potential for a given year,
!      day of year,UT, F10.7, B_z(K_p)
! - low/midlatitudes electric potential is from an empirical model from
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! - the transition zone is smoothed
! - output is the horizontal global electric field in magnetic coordinates direction
!  at every magnetic local time grid point expressed in degrees (0 deg-0MLT; 360deg 24 MLT)
!
! input 
!      integer :: iday,     ! day number of year
!                 iyear     ! year
!      real:: ut,       ! universal time 
!                 F10.7,    ! solar flux       (see ionosphere module)
!                 bz        ! component of IMF (see ionosphere module)
! output
!      real ::               &
!       ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
!       ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
!
! notes:
!
! - !to be done (commented out): input S_a F10.7/ Kp from WACCM and calculate B_z 
!    from these inputs
! - assume regular geomagnetic grid 
! - uses average year 365.24 days/year 30.6001 day/mo s. Weimer
! - get_tilt works only for iyear >= 1900
! - Weimer model 1996, Dan Weimer (not with the updates from B.Emery)
! - fixed parameters: B_z, B_y units nT  CHANGE THIS
!                     F10.7
! - we assume that the reference height is 300km for the emperical potential model
! - as a first approximation the electric field is constant in height
!   WATCH what is the upper boundary condition in WACCM
! - for all the calculation done here we set the reference height to the same 
!   value as in tiegcm (hr=130km)
! - 12/15/03 input value iseasav : replaced by day -> month and day of month
! - 12/15/03 S_aM calculated according to Scherliess draft paper and added
!   S_aM(corrected) = 90*(S_aM+1) to get variation in fig 1 Scherliess draft
!
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

c     use shr_kind_mod,  only: r8 => shr_kind_r8
c     use physconst,     only: pi
c     use abortutils,    only: endrun
c     use cam_logfile,   only: iulog
   
      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: set_fkflfs

      contains
!-----------------------------------------------------------------------
      subroutine set_fkflfs( fk, fl, fs )
!------------------------------------------------------------------
! Purpose:  set f_-k(day) depending on seasonal flag used for empirical model
!     to calculate the electric potential
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/20/03
!-----------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
      USE module_ff ,ONLY:ff
!-----------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------
      real, intent(out) ::  
     &	fk(0:2),  	                ! f_-k(day) 
     &	fl(-2:2), 	                ! f_l(ut)  
     &	fs(2)		                ! f_s(f10.7) 
!------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      integer  :: lp
      real :: ang
      real :: lon_ut

!------------------------------------------------------------------
! f_-k(day) 
! use factors for iseasav == 0 - Scherliess had iseasav as an input parameter
!------------------------------------------------------------------
      lp = iseasav
      if( iseasav == 0 ) then
        ang   = (day + 9.)*dy2rd
        fk(0) = sqr2*cos( 2.*ang )
        fk(1) = sqr2*cos( ang )
        fk(2) = 1.
      else if( iseasav >= 1 .and. iseasav <= 3 ) then
        fk(0) = ft(lp,0)
        fk(1) = ft(lp,1)
        fk(2) = ft(lp,2)
      else if( iseasav == 4 ) then
        fk(0) =0.
        fk(1) =0.
        fk(2) =1.
      end if

!-----------------------------------------------------------------
! f_l(ut) 
!-----------------------------------------------------------------
      lon_ut = 15.*ut        ! 15.*mlt - xmlon + 69. 
      call ff( lon_ut, 2, fl )                                                 
      if( iutav ) then  	! UT-averaging
     
	ang   = fl(0)
        fl(:) = 0.
        fl(0) = ang
	
      end if

!-----------------------------------------------------------------
! f_s(f10.7)  only fs(1) used  	
!-----------------------------------------------------------------
      fs(1) = 1.
!     fs(2) = S_a			  
      fs(2) = f107d			  

      end subroutine set_fkflfs
!-----------------------------------------------------------------------
      end module module_set_fkflfs
