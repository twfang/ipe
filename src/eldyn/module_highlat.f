!20111107: copied originally from tiegcm1.8_dynamo_lres
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
      module module_highlat
!
      PRIVATE
      PUBLIC :: highlat
      
      contains 
!-----------------------------------------------------------------------
      subroutine highlat
      use params_module,ONLY: kmlonp1
      use dynamo_module,ONLY: kmlat0,phihm,potential_model
      use module_sub_heelis,ONLY: heelis
      use module_colath,ONLY:colath
      implicit none
      integer :: i,j
      
!   
! am 10/04 remove weimer part: first test without any potential model
! Dynamo calls Heelis (heelis.F), Weimer (wei01gcm.F), or neither
!   for high latitude electric potential, according to user-provided
!   "model_potential".
! Get high latitude (Heelis or other) colatitudes, NH pfrac, and poten phihm.
!  If Weimer is used, then theta0,phid etc is changed before use in aurora
!   in dynamics.
!
      print *, 'sub-highlat: potential_model=',potential_model
      if (potential_model == 'HEELIS') then
        call heelis
      else  !  potential_model='NONE'
        do j=1,kmlat0
          do i=1,kmlonp1
	    phihm(i,j) = 0.
          enddo ! i=1,kmlonp1
        enddo ! j=1,kmlat0
        call colath
      endif
      
      end subroutine highlat
!-----------------------------------------------------------------------
      end module module_highlat
