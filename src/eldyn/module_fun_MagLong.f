!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!sep2012: efield.f was separated into each routin for SMS compatibility.
!nm20121012: it looks like this routine is not used!!! (thus not included in Makefile) it is called from the subroutine ROTATE which is also not used. this routine must be related to function get_tilt (weimer96) because it shares the same common block: TRANSDAT.
!--------------------------------------------  
      module module_fun_MagLong
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: MagLong

      contains
!-----------------------------------------------------------------------
	FUNCTION MagLong(MLT)
!
!-----------------------------------------------------------------------
! return magnetic longitude in degrees, given Magnetic Local Time
! assuming that TRANS has been called with the date & time to calculate
! the rotation matrices.
!
! btf 11/06/03:
! Call sub adjust instead of referencing it as a function
!-----------------------------------------------------------------------
!
!c       use shr_kind_mod, only: r8 => shr_kind_r8
        USE module_ADJUST, ONLY: ADJUST
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real MagLong
!
!-------------------------------Commons---------------------------------
!
        real cx, st, ct, am
        COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!
!------------------------------Arguments--------------------------------
!
	REAL MLT
!
!---------------------------Local variables-----------------------------
!
	REAL angle, rotangle
!
!-----------------------------------------------------------------------
!
	RotAngle=CX(7)
	angle=MLT*15.-RotAngle-180.
!       MagLong=ADJUST(angle)
        call adjust(angle)
        MagLong = angle
	RETURN
	END FUNCTION MagLong

!-----------------------------------------------------------------------
      end module module_fun_MagLong
