!20111107: copied originally from GIP:  ML__FIELD_LINE_INTEGRALS
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
      MODULE module_calculate_field_line_integrals
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: calculate_field_line_integrals

      CONTAINS
!-----------
      SUBROUTINE calculate_field_line_integrals (NPTS,in,is,ds,midpoint, &
!nm20111107                                         Apex_d1,Apex_d2,&
     & Apex_BE3, &
     & apex_D,apex_d1d1,apex_d1d2,apex_d2d2,apex_BMAG, &
     & ni_oplus_1d,ni_hplus_1d,ti_oplus_1d,no_plus,o2_plus, &
     & n2_plus,n_plus, &
     & tn,o,o2,n2, &
!nm20111107                                         U_zonal,U_merid,U_vert, &
     & Ue, &
     & electron_density_out, &
     & sigma_ped,sigma_hall, & !output
     & sigma_phph_dsi,sigma_lmlm_msi, &
     & sigma_h,sigma_c, &
     & Kdmph_dsi,Kdmlm &
     &,Je,mp,lp )

      USE module_precision
      USE module_input_parameters,ONLY:sw_3DJ,mype
      USE module_eldyn,ONLY: Ed1_90, Ed2_90  !nm20130828, lpconj
      USE module_physical_constants,ONLY:electron_charge_Coulombs, zero
      USE module_IONNEUT_PLAS,ONLY:IONNEUT_PLAS
      USE module_calc_collfreq,ONLY:calc_collfreq
      implicit none

! Input Arguments:
  INTEGER (KIND=int_prec), intent(in)  ::  NPTS
  INTEGER (KIND=int_prec), intent(in)  ::  in    ! North footpoint of a flux tube
  INTEGER (KIND=int_prec), intent(in)  ::  is    ! South footpoint of a flux tube
  INTEGER (KIND=int_prec), intent(in)  ::  midpoint
  REAL (KIND=real_prec), intent(in)  ::  ds(npts)            !ds 1D [???

!nm20111107  REAL (KIND=real_prec), intent(in)  ::  Apex_d1(3,npts) !d1 of reference above (3.8)
!nm20111107  REAL (KIND=real_prec), intent(in)  ::  Apex_d2(3,npts) !d2 of reference above (3.9)

!nm20130830: Be3 is constant along magnetic field lines!
  REAL (KIND=real_prec), intent(in)  ::  Apex_BE3 !B_e3 of reference above (= Bmag/D),[Tesla] (4.13)
  REAL (KIND=real_prec), intent(in)  ::  Apex_D(npts)
  REAL (KIND=real_prec), intent(in)  ::  Apex_d1d1(npts)
  REAL (KIND=real_prec), intent(in)  ::  Apex_d1d2(npts)
  REAL (KIND=real_prec), intent(in)  ::  Apex_d2d2(npts)
  REAL (KIND=real_prec), intent(in)  ::  Apex_BMAG(npts)

  REAL (KIND=real_prec), intent(in)  ::  ni_oplus_1d(npts)  !O+ densities [m-6]
  REAL (KIND=real_prec), intent(in)  ::  ni_hplus_1d(npts)  !H+ densities [m-6]
  REAL (KIND=real_prec), intent(in)  ::  ti_oplus_1d(npts)  !O+ temperatures [K]
  REAL (KIND=real_prec), intent(in)  ::  no_plus(npts)  !NO+ density [m-6]
  REAL (KIND=real_prec), intent(in)  ::  o2_plus(npts)  !O2+ density [m-6]
  REAL (KIND=real_prec), intent(in)  ::  n2_plus(npts)  !N2+ density [m-6]
  REAL (KIND=real_prec), intent(in)  ::  n_plus(npts)  !N+ density [m-6]
  REAL (KIND=real_prec), intent(in)  ::  tn(npts)         !Tn [K]
  REAL (KIND=real_prec), intent(in)  ::  o(npts)          !atomic oxygen density []
  REAL (KIND=real_prec), intent(in)  ::  o2(npts)         !molecular oxygen density []
  REAL (KIND=real_prec), intent(in)  ::  n2(npts)         !molecular nitrogen density []

!nm20111107  REAL (KIND=real_prec), intent(in)  ::  U_zonal(npts)    !neutral wind +ggeast  [m/s]
!nm20111107  REAL (KIND=real_prec), intent(in)  ::  U_merid(npts)    !neutral wind +ggsouth
!nm20111107  REAL (KIND=real_prec), intent(in)  ::  U_vert(npts)     !neutral wind +ggup

! Output Arguments: 

  REAL (KIND=real_prec),    intent(out) ::    sigma_phph_dsi(2)      !(5.13) divided by |sin I_m |
  REAL (KIND=real_prec),    intent(out) ::    sigma_lmlm_msi(2)      !(5.14) multiplied by | sin I_m |

  REAL (KIND=real_prec),    intent(out) ::	    sigma_h(2)      !(5.17)
  REAL (KIND=real_prec),    intent(out) ::	    sigma_c(2)      !(5.18)
  REAL (KIND=real_prec),    intent(out) ::         Kdmph_dsi(2)      !(5.19) divided by |sin I_m |
  REAL (KIND=real_prec),    intent(out) ::	    Kdmlm(2)	  !(5.20) plus or minus ????

  REAL (KIND=real_prec),    INTENT(OUT) :: Je(npts,2)


INTEGER (KIND=int_prec), intent(in)  ::  mp
INTEGER (KIND=int_prec), intent(in)  ::  lp
  REAL (KIND=real_prec), INTENT(IN)  ::  Ue(npts,2)        !Un parallel to e1 [m/s]  (5.6)
!nm20130830:  REAL (KIND=real_prec),INTENT(IN)    ::  Ue2(npts)        !Un parallel to e2 [m/s]
!----------------------------Local variables-----------------------------

  INTEGER (KIND=int_prec) ::  ipts     !
  REAL (KIND=real_prec)    ::  effective_temp(npts)   !???
  REAL (KIND=real_prec)    ::  ion_neut_cf(npts)  ! collision frequency []
  REAL (KIND=real_prec)    ::  ion_mass_amu(npts) ! ion mass in [AMU]
  INTEGER (KIND=int_prec) ::  iout   !???
  REAL (KIND=real_prec)    ::  integral513
  REAL (KIND=real_prec)    ::  integral514
  REAL (KIND=real_prec)    ::  integral517
  REAL (KIND=real_prec)    ::  integral518
  REAL (KIND=real_prec)    ::  integral519
  REAL (KIND=real_prec)    ::  integral520
  INTEGER (KIND=int_prec) ::  ihem  !!South(2) & Northern(1) hemisphere
  INTEGER (KIND=int_prec) ::  istart
  INTEGER (KIND=int_prec) ::  istop
  INTEGER (KIND=int_prec) ::  istep

  REAL (KIND=real_prec)    ::  abs_ds     ! |ds(ipts)|
  REAL (KIND=real_prec)    ::  electron_density  ![m-6]
  REAL (KIND=real_prec)    ::  r_factor
  REAL (KIND=real_prec)    ::  sigma_ped_old  !pedersen conductivity [mho/m]
  REAL (KIND=real_prec)    ::  sigma_hall_old  !pedersen conductivity [mho/m]
  REAL (KIND=real_prec)    ::  sigma_ped(npts)  !pedersen conductivity [mho/m]
  REAL (KIND=real_prec)    ::  sigma_hall(npts) !hall conductivity
  REAL (KIND=real_prec)    ::  electron_density_out(npts) !Ne

  REAL (KIND=real_prec)    ::  qe_fac
  REAL (KIND=real_prec)    ::  o_cm3(npts)         !molecular oxygen density []
  REAL (KIND=real_prec)    ::  o2_cm3(npts)         !molecular oxygen density []
  REAL (KIND=real_prec)    ::  n2_cm3(npts)         !molecular oxygen density []
  REAL (KIND=real_prec)    ::  rnu_o2p(npts)  ! [[o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~n2]n(n2)]/w(o2p)
  REAL (KIND=real_prec)    ::  rnu_op(npts)   ! [[op ~o2]n(o2)+[op ~o]n(o)+[op ~n2]n(n2)]/w(op )
  REAL (KIND=real_prec)    ::  rnu_nop(npts)  ! [[nop~o2]n(o2)+[nop~o]n(o)+[nop~n2]n(n2)]/w(nop)
  REAL (KIND=real_prec)    ::  rnu_ne(npts)   ! electron~neutral
 


! na100504:  sinIm not needed
! real(kind=8)    ::  sinlm    ! sin(lam_m)
! real(kind=8)    ::  clm2     ! cos^2(lam_m)
! real(kind=8)    ::  sinIm    ! sin I_m of reference above (3.7)
! real(kind=8)    ::  abs_sinIm  ! | sin I_m |

!nm20111107  REAL (KIND=real_prec)    ::  electron_charge_Coulombs
!nm20111107  electron_charge_Coulombs=1.6022E-19   ! electronic charge [C]

! when sw_3DJ ON
  REAL (KIND=real_prec)    :: ed11,ed21
  INTEGER (KIND=int_prec) ::  lp0
!---------------------------------------------------------

      sigma_ped=0.0  !JFM this is temporary until Naomi figures our what sigma_ped should be set to

      do ipts = in,is
         effective_temp(ipts) = (tn(ipts)+ti_oplus_1d(ipts))/2.
         IF ( effective_temp(ipts)<tn(ipts) ) effective_temp(ipts) = tn(ipts)
      enddo

      iout = 0
! get ion_mass_amu & ion_neut_cf
      CALL IONNEUT_PLAS(NPTS,o,o2,n2, ni_oplus_1d , no_plus, o2_plus, &
     &               effective_temp,ion_neut_cf,ion_mass_amu, &
     &               in,is,iout)

! calculates collision frequencies see TGCM
  o_cm3 (in:is) = o (in:is)*1.e-6   ! convert from #/m3 to #/cm3
  o2_cm3(in:is) = o2(in:is)*1.e-6
  n2_cm3(in:is) = n2(in:is)*1.e-6
!
  call calc_collfreq(o_cm3,o2_cm3,n2_cm3,effective_temp,tn,apex_Bmag, &
                    rnu_o2p,rnu_op,rnu_nop,rnu_ne,in,is,npts)


!nm20140407: ihem inconsistency to get_efiel90km.f was corrected
!   midpoint = (NPTS+1)/2
! separate South(ihem=2) & Northern(ihem=1) hemisphere
      ihem_loop: do ihem=1,2

! (2) southern hemisphere
         if (ihem==2) then
            istart=is
            istop=midpoint 
            istep=-1
  ! (1) Northern hemisphere
         else if (ihem==1) then
            istart=in
            istop=midpoint - 1
            istep=+1
         endif


         integral513=zero
         integral514=zero
         integral517=zero
         integral518=zero
         integral519=zero
         integral520=zero


! get integral
         ipts_loop1: do  ipts=istart, istop, istep

! get pedersen & hall conductivities

!  electron_density = ni_oplus_1d(ipts)+ni_hplus_1d(ipts)+ no_plus(ipts)+o2_plus(ipts) + n2_plus(ipts)+n_plus(ipts)

!  ....just use the O+, NO+ and O2+ for the Electron density here (original equation above)........
!nm20130905: added all the other ion species for completeness
!           electron_density = ni_oplus_1d(ipts) + no_plus(ipts) + o2_plus(ipts) &
!    &                       + ni_hplus_1d(ipts) + n_plus(ipts)  + n2_plus(ipts)
            electron_density = ni_oplus_1d(ipts) + no_plus(ipts) + o2_plus(ipts)  


            if (electron_density.gt.1.e-10) then


               r_factor         = ion_mass_amu(ipts)*ion_neut_cf(ipts)/electron_charge_Coulombs/apex_BMAG(ipts)
               sigma_ped_old  = electron_density*electron_charge_Coulombs/apex_BMAG(ipts)*r_factor/(1+(r_factor*r_factor))
               sigma_hall_old = sigma_ped(ipts)*r_factor

               qe_fac = electron_charge_Coulombs/apex_BMAG(ipts)

! densities #/m^3
! Pedersen conductivity
        sigma_ped(ipts) = qe_fac* &
    &         ((ni_oplus_1d(ipts)*rnu_op (ipts)/(1.+rnu_op (ipts)**2))+ &
    &          (o2_plus(ipts) *rnu_o2p(ipts)   /(1.+rnu_o2p(ipts)**2))+ &
    &          (no_plus(ipts) *rnu_nop(ipts)   /(1.+rnu_nop(ipts)**2))+ &
    &          (electron_density*rnu_ne (ipts) /(1.+rnu_ne(ipts)**2)))

! Hall conductivity
        sigma_hall(ipts) = qe_fac* &
    &          (electron_density /(1.+rnu_ne (ipts)**2)- &
    &           ni_oplus_1d(ipts)/(1.+rnu_op (ipts)**2)- &
    &           o2_plus(ipts)    /(1.+rnu_o2p(ipts)**2)- &
    &           no_plus(ipts)    /(1.+rnu_nop(ipts)**2))

        electron_density_out(ipts)=electron_density

! get neutral wind vectors
! Ue1=d1*u    :e1:+east (5.6)
!nm20111107  Ue1(ipts)= apex_d1(1,ipts)*U_zonal(ipts) &
!nm20111107  -apex_d1(2,ipts)*U_merid(ipts) &
!nm20111107  +apex_d1(3,ipts)*U_vert(ipts)
!    Ue1(ipts)= -apex_d1(2,ipts)*50.    ! am 071608 test integrals

! Ue2=d2*u    :e2:+down/equatorward
!nm20111107  Ue2(ipts)= apex_d2(1,ipts)*U_zonal(ipts) &
!nm20111107  -apex_d2(2,ipts)*U_merid(ipts) &
!nm20111107  +apex_d2(3,ipts)*U_vert(ipts)
!    Ue2(ipts)= 0.  ! am 071608 test integrals


               abs_ds = ABS(ds(ipts))

! get integrals

!g
!g  The following integrals all come from page 203 and 204 of the paper.  They are numbered
!g  to match the equations in the paper...
!g  The integral parts are calculated here and then some additional factors are applied below.
!g  Note, however, we ignore the |sin I_m| factor wherever it appears since this is applied
!g  later on within the Dynamo solver...
!g

               integral513 = integral513 + sigma_ped(ipts)*apex_d1d1(ipts)*abs_ds/apex_D(ipts)
               
               integral514 = integral514 + sigma_ped(ipts)*apex_d2d2(ipts)*abs_ds/apex_D(ipts)
               
               integral517 = integral517 + sigma_hall(ipts)*abs_ds

               integral518 = integral518 + sigma_ped(ipts)*apex_d1d2(ipts)*abs_ds/apex_D(ipts)
!if(integral518<0.0) then
!!SMS$ignore begin
!  print*,'JFM1', integral518,  sigma_ped(ipts),apex_d1d2(ipts),abs_ds,  apex_D(ipts),  ipts,  ihem,mype
!!SMS$ignore end
!         JFM1 -2.3854740E-05  3.9203514E-06   -2.6396289E-03  2180.000 0.9456919      425      2    0
!         JFM1 -2.4191562E-02  4.8058820E-07   -1.0475066E-02  10876.00 0.6538433      365      2    0
!         JFM1 -2.4609566E-02  0.0000000E+00    1.7838729E-06  4173310. 8.8848356E-05  213      2    0
!This does not happen for mype=1 !!!!!!!!!!!!
!endif

               integral519 = integral519 + (sigma_ped(ipts)*apex_d1d1(ipts)*Ue(ipts,2)/apex_D(ipts) &
     &               + (sigma_hall(ipts)-sigma_ped(ipts)*apex_d1d2(ipts) &
     &               /apex_D(ipts))*Ue(ipts,1))*abs_ds

               integral520 = integral520 + ( (sigma_hall(ipts)+sigma_ped(ipts)*apex_d1d2(ipts) &
     &               /apex_D(ipts) )*Ue(ipts,2) &
     &               - sigma_ped(ipts)*apex_d2d2(ipts)*Ue(ipts,1)/apex_D(ipts) )*abs_ds

!    integral520=integral520 + ( Ue(ipts,1)/apex_D(ipts) )*abs_ds! am 071608 test integrals


               IF ( sw_3DJ==1 ) THEN
!    IF ( ihem==1 ) THEN !NH
!      lp0=lpconj(lp)
!    ELSE IF ( ihem==2 ) THEN !SH
!      lp0=lp
!    END IF
!eq (5.7)
!note: Ed1_90, Ed2_90, Be3 are constant along magnetic field lines!
                  ed11 = Ed1_90(1,lp,mp) + Ue(ipts,2) * Apex_BE3
                  ed21 = Ed2_90(1,lp,mp) - Ue(ipts,1) * Apex_BE3
                  Je(ipts,1) = sigma_ped(ipts)*Apex_d1d1(ipts) * ed11 & 
                       & +         ( sigma_ped(ipts)*Apex_d1d2(ipts) - sigma_hall(ipts)*apex_D(ipts) ) * ed21
                  !eq (5.8)
                  Je(ipts,2) =( sigma_ped(ipts)*Apex_d1d2(ipts) + sigma_hall(ipts)*apex_D(ipts) ) * ed11 &
                       & +            sigma_ped(ipts)*Apex_d2d2(ipts) * ed21

               END IF !( sw_3DJ==1 ) THEN
            
            else             !if (electron_density<1.e-10) then
               print *, 'sub-cal_FLI: !STOP! INVALID Ne', electron_density,ipts,lp,mp
               STOP
            end if !(electron_density.gt.1.e-10) then !(stops electron density = 0.0 producing NaNs)

         enddo  ipts_loop1 !: do  ipts=istart, istop



! inputs to the dynamo solver

!g
!g  Integrals 5.13 and 5.14 do not !include the |sin I_m| or 1/|sin I_m| factors respectively
!g  as these are dealt with in the dynamo module....
!g
         sigma_phph_dsi(ihem) = integral513   !(5.13) divided by |sin I_m |
         sigma_lmlm_msi(ihem) = integral514   !(5.14) multiplied by | sin I_m |
!g
!g  Integrals 5.17 and 5.18.....
!g
         sigma_h(ihem) = integral517       !(5.17)
         sigma_c(ihem) = integral518       !(5.18)
!g
!g  integral 5.19 is multiplied by BE3.  However we have not multiplied by |sin I_m| because
!g  this is done within the dynamo module itself (I've said this enough yeh ?)
!g
         Kdmph_dsi(ihem) =  Apex_BE3*integral519  !(5.19) divided by |sin I_m |
!g
!g  The following is equation 5.20.  The integral is multiplied by BE3.  There is also a minus
!g  sign for the northern hemisphere part (and a plus sign for the southern hemisphere part).
!g  The +- statement is that the upper (lower) sign applies to the northern (southern) magnetic
!g  hemisphere.  This statement is written on page 200 of the paper - just below equation 3.22
!g
!(2) SH
         if(ihem == 2) then
            Kdmlm(ihem) = + Apex_BE3*integral520  !(5.20) plus for southern hemi
!(1) NH
         else if(ihem == 1) then
            Kdmlm(ihem) = - Apex_BE3*integral520  !(5.20) minus for northern hemi
         end if

      enddo ihem_loop  !: do ihem=1,2




END SUBROUTINE calculate_field_line_integrals
END      MODULE module_calculate_field_line_integrals

