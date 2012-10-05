!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_prep_pnm
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
      public :: prep_pnm

      contains

                                                                      



      subroutine prep_pnm
!-----------------------------------------------------------------      
! Purpose: constant factors for normalized associated Legendre polynomial P_n^m
!          Ref.: Richmond J.Atm.Ter.Phys. 1974
!
! Method:
!   PmoPmmo(m) = sqrt(1+1/2m)
!   R_n^m      = sqrt[ (n^2-m^2)/(4n^2-1) ]
!
! Author: A. Maute Nov 2003  am 11/18/03
!-----------------------------------------------------------------     
!nm20121003
      USE efield !,ONLY:
      implicit none                

!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------
      integer  :: mp, m, n
      real :: xms, xns, den

      do mp = 1, mmp            ! m+1 = 1,mm+1                                     
	m = mp - 1                                               
	xms = m*m                                                
	if( mp /= 1 ) then
           pmopmmo(m) = sqrt( 1. + .5/M )
	end if
	do n = m,nm      ! n = m,N                                     
	  xns    = n*n                                       
	  den    = max(4.*xns - 1.,1.)
	  r(n,m) = sqrt( (xns  - xms)/den )
	end do                 
      end do 

      end subroutine prep_pnm                                                                         



      end module module_prep_pnm

