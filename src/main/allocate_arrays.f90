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
!
      SUBROUTINE allocate_arrays ( switch )
      USE module_precision
      USE module_IPE_dimension,ONLY: NMP,NLP,ISTOT
      USE module_FIELD_LINE_GRID_MKS,ONLY: &
     & plasma_grid_3d,plasma_3d &
     &,apexD ,apexE, VEXBup &
     &,Be3, Pvalue, JMIN_IN, JMAX_IS,hrate_cgs_save &
     &,mlon_rad, plasma_grid_Z, plasma_grid_GL
  
      USE module_input_parameters,ONLY: sw_neutral_heating_flip
      IMPLICIT NONE
      INTEGER (KIND=int_prec),INTENT(IN) :: switch
      INTEGER (KIND=int_prec) :: stat_alloc

! (0) ALLOCATE arrays
IF ( switch==0 ) THEN
print *,'ALLOCATing ARRAYS'
      ALLOCATE ( &
!---field line grid
!    &    plasma_grid_3d(NPTS2D, 1:NMP,6) &
!    &    plasma_grid_Z( NPTS2D             ) &
!    &,   plasma_grid_GL(NPTS2D             ) &
!---
!    &,        apexD(3:3,NPTS2D, 1:NMP,3) &
!    &,        apexE(2,NPTS2D,NMP,3) &
!---
     &           Be3(2,1:NMP,NLP) &
     &,       Pvalue(        NLP) &
     &,      JMIN_IN(        NLP) &
     &,      JMAX_IS(        NLP) &
!---
     &,     mlon_rad(  NMP+1    ) &
!---plasma
!    &,    plasma_3d(ISTOT,NPTS2D,1:NMP) &
!dbg20120501     &,    plasma_3d(   1:NMP,NLP) &
!dbg20120501     &,    plasma_3d4n( NPTS2D, 1:NMP) &
     &,    VEXBup(      1:NMP,NLP) &
     &,STAT=stat_alloc         )
 
      IF ( stat_alloc==0 ) THEN
        print *,'ALLOCATion SUCCESSFUL!!!'
      ELSE !stat_alloc/=0
        print *,switch,"!STOP! ALLOCATION FAILD!:",stat_alloc
        STOP
      END IF

! (1) DEALLOCATE arrays
ELSE IF ( switch==1 ) THEN
print *,'DE-ALLOCATing ARRAYS'
! field line grid
      DEALLOCATE ( &
     &    plasma_grid_3d &
!---
     &,        apexD &
!dbg20110923     &,        apexE &
!---
     &,          Be3 &
     &,       Pvalue &
     &,      JMIN_IN &
     &,      JMAX_IS &
!---
     &,      mlon_rad &
!---plasma
     &,  plasma_3d   &
!dbg20120501     &,  plasma_3d4n &
     &,  VEXBup      &
     &,STAT=stat_alloc         )
 
      IF ( stat_alloc==0 ) THEN
        print *,'DE-ALLOCATion SUCCESSFUL!!!'
      ELSE !/=0
        print *, ALLOCATED( plasma_grid_3d )
        print *,switch,"!STOP! DEALLOCATION FAILD!:",stat_alloc
        STOP
      END IF


!---neutral heating
      IF ( sw_neutral_heating_flip==1 ) THEN
         DEALLOCATE ( hrate_cgs_save &
              &,  STAT=stat_alloc         )
         IF ( stat_alloc==0 ) THEN
            print *,'DE-ALLOCATion SUCCESSFUL!!! NHEAT'
         ELSE !/=0
            print *, ALLOCATED( hrate_cgs_save )
            print *,switch,"!STOP! DEALLOCATION FAILD!: NHEAT",stat_alloc
            STOP
         END IF
      END IF !( sw_neutral_heating_flip==1 ) THEN


END IF !( switch==1 ) THEN

      END SUBROUTINE allocate_arrays
