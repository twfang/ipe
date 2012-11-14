!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_read_acoef
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

!     use shr_kind_mod,  only: r8 => shr_kind_r8
!     use physconst,     only: pi
!     use abortutils,    only: endrun
!     use cam_logfile,   only: iulog
   
      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: read_acoef

      contains
!-----------------------------------------------------------------------
      subroutine read_acoef (efield_lflux_file, efield_hflux_file)
!----------------------------------------------------------------     
! Purpose:  
!    1. read in coefficients A_klmn^lf for solar cycle minimum and
!      A_klmn^hf for maximum 
!    2. adjust S_a (f107d) such that if S_a<80 or S_a > 220 it has reasonable numbers
!      S_aM = [atan{(S_a-65)^2/90^2}-a90]/[a180-a90]
!      a90  = atan [(90-65)/90]^2
!      a180 = atan [(180-65)/90]^2
!    3. inter/extrapolation of the coefficient to the actual flux which is
!      given by the user
!      A_klmn = S_aM [A_klmn^hf-A_klmn^lf]/90. + 2*A_klmn^lf-A_klmn^hf
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/19/03
!---------------------------------------------------------------

!     use ioFileMod,     only : getfil
!     use units,         only : getunit, freeunit
!nm20121003
      USE efield !,ONLY:

      character(len=*), intent(in) :: efield_lflux_file
      character(len=*), intent(in) :: efield_hflux_file

      integer  :: i,ios,unit,istat
      character (len=256):: locfn

!SMS$SERIAL(<a_lf,OUT> : default=ignore) BEGIN
!------------------------------------------------------------------    
!  get coefficients file for solar minimum: 
!-----------------------------------------------------------------                                                                   
!     unit     = getunit()
      unit     = 11
!     call getfil( efield_lflux_file, locfn, 0 )
      locfn=efield_lflux_file

!------------------------------------------------------------------    
! open datafile with coefficients A_klnm
!------------------------------------------------------------------     
!     write(iulog,*) 'read_acoef: open file ',trim(locfn),
!    &' unit ',unit
!dbg20110906
      print *,'unit',unit
      open(unit=unit,file=trim(locfn),status = 'old',iostat = ios)
!     if(ios.gt.0) then
!     write(iulog,*) 
!    &'read_acoef: error in opening coeff_lf file',
!    &' unit ',unit
!       call endrun
!     end if

!----------------------------------------------------------------------------                                                                   
! read datafile with coefficients A_klnm
!--------------------------------------------------------------------   
!     write(iulog,*) 'read_acoef: read file ',trim(locfn),' unit ',unit
      read(unit,*,iostat = ios) a_lf
!     if(ios.gt.0) then
!     write(iulog,*) 
!    &'read_acoef: error in reading coeff_lf file',' unit ',unit
!       call endrun
!     end if

!--------------------------------------------------------------------  
! close & free unit      
!--------------------------------------------------------------------  
      close(unit)
!     call freeunit(unit)
!     write(iulog,*) 'read_acoef: free unit ',unit
!SMS$SERIAL END

!SMS$SERIAL(<a_hf,OUT> : default=ignore) BEGIN

!--------------------------------------------------------------------  
!  get coefficients file for solar maximum: 
!--------------------------------------------------------------------
!     unit     = getunit()
      unit     = 10
!     call getfil( efield_hflux_file, locfn, 0 )
      locfn= efield_hflux_file

!-------------------------------------------------------------------
! open datafile with coefficients A_klnm
!------------------------------------------------------------------
!     write(iulog,*) 'read_acoef: open file ',trim(locfn),' unit ',unit
      open(unit=unit,file=trim(locfn),status = 'old',iostat = ios)
!     if(ios.gt.0) then
!      write(iulog,*) 
!    &'read_acoef: error in opening coeff_hf file',' unit ',unit
!       call endrun
!     end if

!-----------------------------------------------------------------
! read datafile with coefficients A_klnm
!----------------------------------------------------------------
!     write(iulog,*) 'read_acoef: read file ',trim(locfn)
      read(unit,*,iostat = ios) a_hf
!     if(ios.gt.0) then
!      write(iulog,*) 
!    &'read_acoef: error in reading coeff_hf file',' unit ',unit
!       call endrun
!     end if

!---------------------------------------------------------------
! close & free unit      
!-------------------------------------------------------------- 
      close(unit)
!     call freeunit(unit)
!     write(iulog,*) 'read_acoef: free unit ',unit
!SMS$SERIAL END

      end subroutine read_acoef
!-----------------------------------------------------------------------
      end module module_read_acoef
