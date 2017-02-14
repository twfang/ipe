!20150717: copied originally from GIP: calc_collfreq
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
      MODULE module_calc_collfreq
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: calc_collfreq

      CONTAINS
!--------------------
      SUBROUTINE calc_collfreq(o1_cm3,o2_cm3,n2_cm3,tnti,tn,apex_Bmag,  &
                    rnu_o2p,rnu_op,rnu_nop,rnu_ne,in,is,npts)
      USE module_precision
      IMPLICIT NONE
      integer(kind=int_prec),intent(in):: in,is,NPTS
      real(kind=real_prec), intent(in)::o1_cm3(npts),o2_cm3(npts),      &
     &                                  n2_cm3(npts),                   &
     &         tnti(npts),tn(npts),apex_Bmag(npts)
      real(kind=real_prec), intent(out):: rnu_o2p(npts),rnu_op(npts),   &
     &         rnu_nop(npts),rnu_ne(npts)
!
      real(kind=real_prec),parameter :: colfac = 1.5,                   &
     &      qeoNao10 =  9.6489E7,                                       &        ! qe/N_a*1000 [C/mol]
     &      qeomeo10 =  1.7588028E11        ! qe/m_e*1000 [C/g]

      real(kind=real_prec),parameter :: rmass_nop = 30.! [g/mol]
      real(kind=real_prec),parameter :: rmassinv_nop=1./rmass_nop ! inverted rmass
      real(kind=real_prec),parameter :: rmass_o1 = 16.! [g/mol]
      real(kind=real_prec),parameter :: rmassinv_o1=1./rmass_o1 ! inverted rmass
      real(kind=real_prec),parameter :: rmass_o2 = 32.! [g/mol]
      real(kind=real_prec),parameter :: rmassinv_o2=1./rmass_o2 ! inverted rmass

! local
      integer(kind=int_prec):: ipts
      real(kind=real_prec) :: te, sqrt_te,                              &
     &      omega_op, omega_o2p, omega_nop,                             &
     &      omega_op_inv, omega_o2p_inv,                                &
     &      omega_nop_inv, omega_e, omega_e_inv

      real(kind=8),dimension(npts)::                                    &
           rnu_o2p_o2, & ! O2+ ~ O2 collision freq (resonant, temperature dependent)
           rnu_op_o2 , & ! O+  ~ O2 collision freq (non-resonant)
           rnu_nop_o2, & ! NO+ ~ O2 collision freq (non-resonant)
!
           rnu_o2p_o, &  ! O2+ ~ O  collision freq (non-resonant)
           rnu_op_o , &  ! O+  ~ O  collision freq (resonant, temperature dependent)
           rnu_nop_o, &  ! NO+ ~ O  collision freq (non-resonant)
!
           rnu_o2p_n2, & ! O2+ ~ N2 collision freq (non-resonant)
           rnu_op_n2 , & ! O+  ~ N2 collision freq (non-resonant)
           rnu_nop_n2    ! NO+ ~ N2 collision freq (non-resonant)

! Ion-neutral momentum transfer collision frequencies [cm^3/s]: check
! units!!!!!
      do  ipts = in , is
!
! gyrofrequencies: omega_i = eB/m_i  [1/s]
!                  omega_e = eB/m_e  [1/s]
! with qeoNao10 = e/Na [C/mol g/kg ]
!      qeomeo10 = e/m_e [C/g g/kg ]
! 1000 in qeoNao10 and qeomeo10 for conversion from 1/g to 1/kg
!
            omega_op     = qeoNao10*apex_Bmag(ipts)*rmassinv_o1
            omega_o2p    = qeoNao10*apex_Bmag(ipts)*rmassinv_o2
            omega_nop    = qeoNao10*apex_Bmag(ipts)*rmassinv_nop
            omega_op_inv = 1./omega_op
            omega_o2p_inv= 1./omega_o2p
            omega_nop_inv= 1./omega_nop
            omega_e      = qeomeo10*apex_Bmag(ipts)
            omega_e_inv  = 1./omega_e
!
! O2 collision frequencies:
           rnu_o2p_o2(ipts) = 2.59E-11*sqrt(tnti(ipts))*  &! O2+ ~ O2 (resonant)
            (1.-0.073*log10(tnti(ipts)))**2
           rnu_op_o2 (ipts) = 6.64E-10                  ! O+  ~ O2
           rnu_nop_o2(ipts) = 4.27E-10                  ! NO+ ~ O2
!
! O collision frequencies:
           rnu_o2p_o(ipts) = 2.31E-10                   ! O2+ ~ O
           rnu_op_o (ipts) = 3.67e-11*sqrt(tnti(ipts))* &  ! O+  ~ O (resonant)
            (1.-0.064*log10(tnti(ipts)))**2*colfac
           rnu_nop_o(ipts) = 2.44E-10                   ! NO+ ~ O
!
! N2 collision frequencies:
           rnu_o2p_n2(ipts) = 4.13E-10                  ! O2+ ~ N2
           rnu_op_n2 (ipts) = 6.82E-10                  ! O+  ~ N2
           rnu_nop_n2(ipts) = 4.34E-10                  ! NO+ ~ N2

! collision frequency nu_in for each ion [1/s]
!    by multiplying with neutral number density [1/cm^3] and sum over
!    neutrals
! nu_in is divided by gyrofrequency omega_i
! nu_in/omega_i [-]:
! rnu_o2p = [[o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~n2]n(n2)]/w(o2p)
! rnu_op  = [[op ~o2]n(o2)+[op ~o]n(o)+[op ~n2]n(n2)]/w(op )
! rnu_nop = [[nop~o2]n(o2)+[nop~o]n(o)+[nop~n2]n(n2)]/w(nop)
           rnu_o2p(ipts) = (rnu_o2p_o2(ipts)*o2_cm3(ipts) + &
                           rnu_o2p_o (ipts)*o1_cm3(ipts) + &
                           rnu_o2p_n2(ipts)*n2_cm3(ipts))*omega_o2p_inv
           rnu_op (ipts) = (rnu_op_o2 (ipts)*o2_cm3(ipts) + &
                           rnu_op_o  (ipts)*o1_cm3(ipts) + &
                           rnu_op_n2 (ipts)*n2_cm3(ipts))*omega_op_inv
           rnu_nop(ipts) = (rnu_nop_o2(ipts)*o2_cm3(ipts) + &
                           rnu_nop_o (ipts)*o1_cm3(ipts) + &
                           rnu_nop_n2(ipts)*n2_cm3(ipts))*omega_nop_inv
!
! neutral~electron collision frequency (from Banks & Kockards) nu_en
! divided by gyrofrequency omega_2:
! nu_en/omega_e [-]
!
            te = tn(ipts)       ! Te is not available approxilate Te = Tn
            sqrt_te = sqrt(te)

            rnu_ne(ipts) =  &
             (2.33e-11*n2_cm3(ipts)*te     *(1.-1.21e-4*te     )+ &
              1.82e-10*o2_cm3(ipts)*sqrt_te*(1.+3.60e-2*sqrt_te)+ &
              8.90e-11*o1_cm3(ipts)*sqrt_te*(1.+5.70e-4*te    ))* &
              omega_e_inv
!
! 6/2/06 btf: Multiply rnu_ne by 4, as per Richmond:
!
! The effective electron-neutral collision frequency is increased in
! an an hoc manner by a factor of 4 in order for the model to produce
! electric fields and currents below 105 km that agree better with
! observations, as recommended by Gagnepain et al. (J. Atmos. Terr.
! Phys., 39, 1119-1124, 1977).
!
           rnu_ne(ipts) = rnu_ne(ipts)*4.
  
      enddo ! loop over points in, is


      END SUBROUTINE calc_collfreq
      END MODULE module_calc_collfreq
