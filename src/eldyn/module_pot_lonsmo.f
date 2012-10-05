!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_pot_lonsmo
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
      public :: pot_lonsmo

      contains
!-----------------------------------------------------------------------
      subroutine pot_lonsmo( pot, idlon ) 
!-------------------------------------------------------------------
! Purpose: smoothing in longitude of potential
!
! Method:  weighted smoothing in longitude
!          assume regular grid spacing
!
! Author: A. Maute Nov 2003  am 11/20/03
!-------------------------------------------------------------------
!nm20121003
      USE efield !,ONLY:
!-------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      integer, intent(in)     :: idlon
      real, intent(inout) :: pot(0:nmlon,0:nmlat)

!-------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
      integer  :: ilon, ilat, id, iabs
      real :: wgt, del
      real :: w(-idlon:idlon)
      real :: pot_smo(0:nmlath) ! temp array for smooth. potential
      real :: tmp(-idlon:nmlon+idlon) ! temp array for smooth. potential

!-------------------------------------------------------------------
! weighting factors (regular grid spacing) 
!-------------------------------------------------------------------
      wgt = 0.
      do id = -idlon,idlon
        del   = abs(id)*dlonm	! delta lon_m
        w(id) = 1./(del + 1.)
        wgt   = wgt + w(id)
      end do
        wgt = 1./wgt

!-------------------------------------------------------------------
! averaging     
!-------------------------------------------------------------------
!     do ilon = 0,nmlon
!       do ilat = 0,nmlath
!         pot_smo(ilat) = 0.
!         do id = -idlon,idlon	                  !  org. was degree now grid points
!           iabs = ilon + id
!           if( iabs > nmlon ) then
!              iabs = iabs - nmlon ! test if wrap around
!           end if
!           if( iabs < 0 ) then
!              iabs = iabs + nmlon ! test if wrap around
!           end if
!           pot_smo(ilat) = pot_smo(ilat) + w(id)*pot(iabs,ilat)
!         end do
!         pot_smo(ilat)  = pot_smo(ilat)*wgt
!         pot(ilon,ilat) = pot_smo(ilat)       ! copy back into pot 
!         pot(ilon,nmlat-ilat) = pot_smo(ilat) ! copy back into pot    
!       end do   
!     end do

!$omp parallel do private(ilat,ilon,tmp)
      do ilat = 0,nmlath
          tmp(0:nmlon)             = pot(0:nmlon,ilat)
          tmp(-idlon:-1)           = pot(nmlon-idlon:nmlon-1,ilat)
          tmp(nmlon+1:nmlon+idlon) = pot(1:idlon,ilat)
          do ilon = 0,nmlon
      pot(ilon,ilat) = dot_product( tmp(ilon-idlon:ilon+idlon),w )*wgt
      pot(ilon,nmlat-ilat) = pot(ilon,ilat)
          end do   
      end do
      
      end subroutine pot_lonsmo
!-----------------------------------------------------------------------
      end module module_pot_lonsmo
