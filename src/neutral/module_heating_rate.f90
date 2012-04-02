!nm20111118: heating rate output available!!!
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
MODULE module_heating_rate
      USE module_precision
      USE module_IPE_dimension,ONLY: NPTS2D,NMP0,NMP1
      IMPLICIT NONE
!
      !.. EHT_cgs(3,J) = e heating rate, EHT_cgs(1,J) = ion heating rate, EHT_cgs(2,J) unused (eV cm-3 s-1)
!20110516: temporary moved to module plasma for the new structure for mpi test
!      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP), PUBLIC :: EHT_cgs
      !.. TE_TI_k(3,J) = Te, TE_TI_k(2,J) = Ti = TE_TI_k(2,J) [kelvin]
!dbg20110927      REAL(KIND=real_prec), DIMENSION(NPTS2D, NMP0:NMP1), PUBLIC :: NHEAT_cgs !.. Neutral heating rate (eV/cm^3/s) 
!tmp20110404: neutral heating rate in eV kg-1 s-1
!t      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP), PUBLIC :: NHEAT_mks !.. Neutral heating rate (eV/kg/s) 
      REAL(KIND=real_prec), DIMENSION(:,:),ALLOCATABLE,PUBLIC :: hrate_cgs_save !.. each component of the Neutral heating rate (eV/cm^3/s) DIM(7,NPTS2D,NMP0:NMP1)
!nm20111118: moved to module_FIELD_LINE_GRID_MKS
!nm20111118:      INTEGER,PUBLIC :: mp_save,lp_save !mp,lp values to be referred to from outside of subroutine plasma


      PRIVATE
      PUBLIC :: get_neutral_heating_rate

      CONTAINS
!-----------------------------
! convert the neutral heating rate from eV/cm3/s --> eV/kg/s -->J/kg/s
      SUBROUTINE get_neutral_heating_rate ( hrate_mks )
      USE module_NEUTRAL_MKS,ONLY: ON_m3,HN_m3,N2N_m3,O2N_m3,HE_m3,N4S_m3
      USE module_physical_constants,ONLY:mass_kg,AMU
      USE module_unit_conversion,ONLY:M3_to_CM3,eV2J
      USE module_FIELD_LINE_GRID_MKS,ONLY: mp_save,plasma_grid_Z
      IMPLICIT NONE
!------------------------
      REAL(KIND=real_prec), DIMENSION(7,NPTS2D),INTENT(OUT) :: hrate_mks !.. each component of the Neutral heating rate (eV/kg/s) 
      REAL(KIND=real_prec) :: total_rho
      INTEGER(KIND=int_prec) :: i,jth

      i_loop: DO i=1,NPTS2D
        IF(plasma_grid_Z(i)>=80.E+03.AND.plasma_grid_Z(i)<=700.E+03) THEN
! rho[kg cm-3] = n[cm-3] * mass[kg]
         total_rho = ( &
              &     ON_m3(i,mp_save) * mass_kg(1) &
              &  +  HN_m3(i,mp_save) * mass_kg(2) &
              &  + N2N_m3(i,mp_save) * mass_kg(3) &
              &  + O2N_m3(i,mp_save) * mass_kg(4) &
              &  +  HE_m3(i,mp_save) * mass_kg(5) &
              &  + N4S_m3(i,mp_save) * mass_kg(6) &
              & ) * M3_to_CM3 * AMU

! convert unit: eV/cm3/s --> J/cm3/s --> J/kg/s
 !t     NHEAT_mks(1:NPTS2D,NMP0:NMP1) = NHEAT_cgs(1:NPTS2D,NMP0:NMP1) * eV2J / total_rho(1:NPTS2D,NMP0:NMP1) 
! each component of the heating rate
          IF ( total_rho>0.0 ) THEN
            hrate_mks(1:7,i) = hrate_cgs_save(1:7,i) * eV2J / total_rho
          ELSE

print *,'!STOP! INVALID total_rho',total_rho,i,mp_save
STOP
!  hrate_mks(1:7,i) = 0.0

          END IF ! ( total_rho>0.0 ) THEN

        ELSE    ! IF(plasma_grid_Z(i)<80.or.plasma_grid_Z(i)>700) THEN

          hrate_mks(1:7,i) = 0.0

        END IF  !(plasma_grid_Z(i)>=80.AND.plasma_grid_Z(i)<=700) THEN      

      END DO i_loop ! i=1,NPTS2D
    END SUBROUTINE get_neutral_heating_rate
  
END MODULE module_heating_rate
