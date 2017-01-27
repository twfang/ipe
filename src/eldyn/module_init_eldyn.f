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
      USE module_eldyn,only :j0,j1,Ed1_90,Ed2_90,coslam_m               &
     &, plas_fli !t,Je_3d
      USE module_input_parameters,only: sw_eldyn
      use module_init_cons,only:init_cons 
      use module_init_heelis,only:init_heelis 
      use module_nc_create,only:nc_create
!t      use module_readin_netcdf,only:readin_netcdf
      use module_readin_ascii,only:readin_ascii
      use read_module,only:input_type
      IMPLICIT NONE
      integer :: status
!20120304:      CHARACTER(len=*),PARAMETER :: path='~/sandbox/efield/'

      print *,'begin init_eldyn: sw_eldyn=', sw_eldyn

!1: self-consistent electrodynamic solver
!t      IF ( sw_eldyn==0 ) THEN 

          print *,'sub-init_eldyn: calling init_cons'
          CALL init_cons

          print *,'sub-init_eldyn: calling init_heelis'
          CALL init_heelis

          print *,'sub-init_eldyn: calling nc_create'
          CALL nc_create

          ! read in integrals                                                             
          if(input_type == 'NETCDF') then
!            print *,'sub-init_eldyn: calling readin_netcdf'
!            call readin_netcdf
          elseif(input_type == 'ASCII') then
            print *,'sub-init_eldyn: calling readin_ascii'
            call readin_ascii
          else
            write(6,*)'Did not recognize input_type=',input_type
            stop 'couple'
          endif



  
!2: WACCM empirical electric field model
!t      ELSE IF ( sw_eldyn==1 ) THEN 
        print *,'allocating Ed1_90etc'

      allocate( j0      (2,NLP    ),                                    &
     &          j1      (2,NLP    ),                                    &
     &          coslam_m(2,NLP    ),                                    &
     &          Ed1_90  (2,NLP,NMP),                                    &
     &          Ed2_90  (2,NLP,NMP),                                    &
!nm20131025:...ideally, the array should be located for sw_eldyn=0
     &          plas_fli (2,NLP,NMP,6),                                 &
!t     &          sigma_ped(     2,NLP,NMP),                              &
!t     &          sigma_hall(    2,NLP,NMP),                              &
     &          STAT=status       )
      if(status /=0) then
        print*,'Allocation failed in module_init_eldyn',status
        print*,'Stopping in module_init_eldyn'
        stop
      endif
      CALL efield_init( 'coeff_lflux.dat',                              &
     &                  'coeff_hflux.dat',                              &
     &                  'wei96.cofcnts'   )

!t      END IF !( sw_eldyn==0 ) THEN 

      print *,'END sub-init_eld'
      END SUBROUTINE init_eldyn
!---
      END MODULE module_init_eldyn
