!t #include "defs.h"
!
      module cons_module
      use params_module,only: kmlat,kmlon,kmlonp1
       
      implicit none

!pi: set with 4*atan(1)    C(110)
!dtr: degrees-to-radians (pi/180.) 
!rtd: radians-to-degrees (180./pi)    
      real ::                                                           &
     &  pi,                                                             &
     &  dtr,                                                            &
     &  rtd
      real,parameter ::                                                 &
     &  pi_dyn=3.14159265358979312 ! pi used in dynamo calculations
!re: earth radius (cm)  
!h0 =9.0e6, r0 =re+h0: use mean earth radius              
      real,parameter ::                                                 &
     &  re   = 6.37122e8,                                               &
     &  h0 =9.0e6, r0 =re+h0
!dlat,dlonm: grid spacing
!xlatm: magnetic latitudes (radians)
!xlonm:! magnetic longitudes (radians)
!xlatm_deg:! magnetic latitudes (degree)
!xlonm_deg:! magnetic longitudes (degree)
!rcos0s(kmlat):! cos(theta0)/cos(thetas)
!dt1dts(kmlat):  ! dt0dts/abs(sinim) (non-zero at equator)
      real, public ::                                                   &
     &  dlatm,                                                          &
     &  dlonm,                                                          &
     &  xlatm(kmlat),                                                   &
     &  xlonm(kmlonp1),                                                 &
     &  xlatm_deg(kmlat),                                               &
     &  xlonm_deg(kmlonp1),                                             &
     &  rcos0s(kmlat),                                                  &
     &  dt1dts(kmlat)
!
! Critical colatitude limits for use of Heelis potential in dynamo:
      real,parameter ::                                                 &
     &  crit(2) = (/0.523598775, 0.61086524/)  ! plasmasphere has zero  
!     &  crit(2) = (/0.261799387, 0.523598775/) ! original values
!nm031407:     &  crit(2) = (/0.523598775, 0.525298775/)  !nm041106: test
                                          ! conductances aboce &lam&>60 deg therefore I set the pcb to 60deg
                                          ! and the auroral zone equator boundary to 55 deg
	! now above &lam&>60 set the conductances to some value which shouldn't
	! matter and below it's an interpolation between high and low latitude 
	! potential
!nm040607:     &  crit(2) = (/0.5314827561378479, 0.61086524/)  !crit(1) corresponds to max lat gip(83)=59.54828119329674  !nm040607:
	
      integer, parameter :: jyr =1997,                                  &
     &   jday =97
      real, parameter :: jsecs =150300.! JSECS

      integer, public    :: idyn_save(kmlat) !correspondance between lp_plas & lp_dyn
      INTEGER, parameter :: lp_dyn_eq=47     !the lowest latitude index for FLI

!-----------------------------------------------------------------------
      end module cons_module
