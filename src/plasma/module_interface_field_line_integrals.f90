!20111107: copied originally from GIP:  ML__FIELD_LINE_INTEGRALS, ML__IONNEUT_PLAS
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
      MODULE module_interface_field_line_integrals
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: interface_field_line_integrals
      CONTAINS
!-------------
      SUBROUTINE interface_field_line_integrals ( lp_plas,mp,utime,  &
     & sigma_ped_3d,sigma_hall_3d,Ue1_3d,Ue2_3d,Ne_3d )
      USE module_precision
!nm20130906      USE module_IPE_dimension,ONLY: IPDIM
      USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_grid_3d,plasma_grid_GL, apexD, JMIN_IN,JMAX_IS,Be3, apexDscalar,MaxFluxTube, ISL,IBM,east,north,up &
     &, TN_k,ON_m3,O2N_m3,N2N_m3,Un_ms1 &
     &, plasma_3d
      USE module_input_parameters,ONLY:sw_3DJ,mpstop,mype
      USE module_physical_constants,ONLY:rtd
      USE module_IPE_dimension,ONLY: NPTS2D,NMP,NLP
      USE module_calculate_field_line_integrals,ONLY:calculate_field_line_integrals
      USE cons_module,ONLY:xlatm,idyn_save 
!     USE module_PLASMA,ONLY:sigma_ped_3d,sigma_hall_3d,Ue1_3d,Ue2_3d

!nm20140731      USE module_eldyn,ONLY:lp_dyn_min
!nm20150330      USE module_calculate_ylatm1,ONLY:calculate_ylatm1
      USE module_save2fli_array,ONLY:save2fli_array
!t      USE params_module,ONLY:nmlat1
      IMPLICIT NONE
!---argument
INTEGER (KIND=int_prec), intent(in)  ::  mp
INTEGER (KIND=int_prec), intent(in)  ::  lp_plas
INTEGER (KIND=int_prec), INTENT(IN)  ::  utime !universal time [sec]
!---local
INTEGER (KIND=int_prec)  ::  in,in1d    ! North footpoint of a flux tube
INTEGER (KIND=int_prec)  ::  is,is1d    ! South footpoint of a flux tube
INTEGER (KIND=int_prec)  ::  midpoint,midpoint1d
REAL (KIND=real_prec)    ::  ds(MaxFluxTube)            !ds 1D [???

!nm20130830: Be3 is constant along a field line!
  REAL (KIND=real_prec)  ::  Apex_BE3 !B_e3 of reference above (= Bmag/D),[Tesla] (4.13)
  REAL (KIND=real_prec)  ::  Apex_D(MaxFluxTube)
  REAL (KIND=real_prec)  ::  Apex_d1d1(MaxFluxTube)
  REAL (KIND=real_prec)  ::  Apex_d1d2(MaxFluxTube)
  REAL (KIND=real_prec)  ::  Apex_d2d2(MaxFluxTube)
  REAL (KIND=real_prec)  ::  Apex_BMAG(MaxFluxTube)

  REAL (KIND=real_prec)  ::  ni_oplus_1d(MaxFluxTube)  !O+ densities [m-6]
  REAL (KIND=real_prec)  ::  ni_hplus_1d(MaxFluxTube)  !H+ densities [m-6]
  REAL (KIND=real_prec)  ::  ti_oplus_1d(MaxFluxTube)  !O+ temperatures [K]
  REAL (KIND=real_prec)  ::  no_plus(MaxFluxTube)  !NO+ density [m-6]
  REAL (KIND=real_prec)  ::  o2_plus(MaxFluxTube)  !O2+ density [m-6]
  REAL (KIND=real_prec)  ::  n2_plus(MaxFluxTube)  !N2+ density [m-6]
  REAL (KIND=real_prec)  ::  n_plus(MaxFluxTube)  !N+ density [m-6]
  REAL (KIND=real_prec)  ::  tn(MaxFluxTube)         !Tn [K]
  REAL (KIND=real_prec)  ::  o(MaxFluxTube)          !atomic oxygen density []
  REAL (KIND=real_prec)  ::  o2(MaxFluxTube)         !molecular oxygen density []
  REAL (KIND=real_prec)  ::  n2(MaxFluxTube)         !molecular nitrogen density []

  REAL (KIND=real_prec)  ::  Ue(MaxFluxTube,2)        !Un parallel to e1(1)/e2(2) [m/s]  (5.6)
!nm20130830:  REAL (KIND=real_prec)  ::  Ue2(MaxFluxTube)        !Un parallel to e2 [m/s]

! Output Arguments: 

  REAL (KIND=real_prec) ::  sigma_ped_1d(MaxFluxTube)  !pedersen conductivity [mho/m]
  REAL (KIND=real_prec) ::  sigma_hall_1d(MaxFluxTube) !hall conductivity
  REAL (KIND=real_prec) ::  sigma_ped_3d(MaxFluxTube,47,NMP)  !pedersen conductivity [mho/m]
  REAL (KIND=real_prec) ::  sigma_hall_3d(MaxFluxTube,47,NMP)  !hall conductivity [mho/m]
  REAL (KIND=real_prec) ::  Ue1_3d(MaxFluxTube,47,NMP) 
  REAL (KIND=real_prec) ::  Ue2_3d(MaxFluxTube,47,NMP) 
  REAL (KIND=real_prec) ::  electron_density_out(MaxFluxTube) !hall conductivity
  REAL (KIND=real_prec) ::  Ne_3d(MaxFluxTube,47,NMP) 

  REAL (KIND=real_prec) ::    sigma_phph_dsi_1d(2)      !(5.13) divided by |sin I_m |
  REAL (KIND=real_prec) ::    sigma_lmlm_msi_1d(2)      !(5.14) multiplied by | sin I_m |
  REAL (KIND=real_prec) ::	    sigma_h_1d(2)      !(5.17)
  REAL (KIND=real_prec) ::	    sigma_c_1d(2)      !(5.18)
  REAL (KIND=real_prec) ::         Kdmph_dsi_1d(2)      !(5.19) divided by |sin I_m |
  REAL (KIND=real_prec) ::	    Kdmlm_1d(2)	  !(5.20) plus or minus ????
! when sw_3DJ=1
  REAL (KIND=real_prec) ::       Je_1d(MaxFluxTube,2) !1D slice of the 3DJ(mp,lp)
!------other local variables----------------------
  INTEGER (KIND=int_prec) ::  CTIPDIM
  INTEGER (KIND=int_prec) ::  NPTS
  INTEGER (KIND=int_prec) ::  i,i1d,jth,jth1,jth2,i_set
  REAL (KIND=real_prec) :: dotprod,d1xd2(3)
!---
  REAL (KIND=real_prec) ::  sigma_ped(NPTS2D,NMP)  !pedersen conductivity [mho/m]
  REAL (KIND=real_prec) ::  sigma_hall(NPTS2D,NMP) !hall conductivity
  REAL (KIND=real_prec) ::  mlat_dyn, mlat_plas
  INTEGER (KIND=int_prec) ::  lp_dyn, imlat_dyn, imlat_plas
  INTEGER (KIND=int_prec),parameter ::  lp_dyn_eq=47 !the lowest latitude index for FLI
!---
!nm20130830: lp-->lp_plas
! if mlat_lp_plas[deg] == mlat_lp_dyn then continue, otherwise, RETURN
      mlat_plas  = 90. - plasma_grid_GL(JMAX_IS(lp_plas),lp_plas)*rtd
      imlat_plas = INT(mlat_plas*10.)
print *,' mp=',mp,' lp_plas=',lp_plas,' mlat_plas=', mlat_plas,' imlat_plas=', imlat_plas

!nm20150330      if ( mp==1.AND.lp_plas==1 ) then !utime==start_time???
!nm20150330         call calculate_ylatm1 ( ) 
if (mp==1.and.lp_plas==1)   print *,'xlatm[deg]=',xlatm*rtd
! STOP
!nm20150330      endif

      lp_dyn_loop: DO lp_dyn=1,lp_dyn_eq !from SH toward eq
         mlat_dyn  = xlatm(lp_dyn)*rtd  ![deg]
         imlat_dyn = INT(mlat_dyn*10.)

         print *,'lp_dyn=',lp_dyn,' mlat_dyn=',mlat_dyn,' imlat_dyn=',imlat_dyn
         !dbg20151107: make sure fli is calculated at lp_plas=NLP
         IF ( lp_plas < NLP ) THEN
            IF ( imlat_dyn > imlat_plas ) THEN
               print *, '(1)FLI not calculated: mp=',mp,' lp_plas=', lp_plas,imlat_plas,imlat_dyn
               EXIT lp_dyn_loop
            ELSE IF ( imlat_dyn < imlat_plas ) THEN
               !d print *,'(2)CYCLE! mlat_dyn=',mlat_dyn,' mlat_plas=', mlat_plas,' imlat_plas=', imlat_plas
               CYCLE lp_dyn_loop
            END IF !( imlat_dyn > imlat_plas ) THEN  !dbg20151107

            print *,'(3) start interface calculating FLI: lp_plas=',lp_plas,lp_dyn,' mlat_dyn=',mlat_dyn
            idyn_save(lp_dyn)=lp_plas  !correspondance between lp_plas & lp_dyn
         ELSE if ( lp_plas == nlp ) then 

print *, '!dbg20151107 make sure fli is calculated at lp_plas=170'
           if ( lp_dyn < lp_dyn_eq )  CYCLE lp_dyn_loop
           print *, lp_dyn, lp_plas
           idyn_save(lp_dyn)=lp_plas  !correspondance between lp_plas & lp_dyn
         END IF !( lp_plas < NLP ) THEN

         IN = JMIN_IN(lp_plas)
         IS = JMAX_IS(lp_plas)
         midpoint = IN + ( IS - IN )/2
         CTIPDIM = IS - IN + 1
         i_loop: DO i=in,is
            i1d=i-in+1
            IF ( i==is ) THEN
               ds(i1d) = plasma_grid_3d(IS,lp_plas,mp,ISL) - plasma_grid_3d(IS-1,lp_plas,mp,ISL)
            ELSE
               ds(i1d) = plasma_grid_3d(i+1,lp_plas,mp,ISL) - plasma_grid_3d(i,lp_plas,mp,ISL)
            END IF

!nm20111107  Apex_d2(1,i1d)=apexD(2,i,mp)%east
!nm20111107  Apex_d2(2,i1d)=apexD(2,i,mp)%north
!nm20111107  Apex_d2(3,i1d)=apexD(2,i,mp)%up

!DOT_PRODUCT( D3(1:3,i,mp), Vn_ms1(1:3,i) ) / SQRT(  DOT_PRODUCT( D3(1:3,i,mp), D3(1:3,i,mp) )  )     
! Un(1)=Ue1=d1*U: positive east, Eq(5.6) 
! Un(2)=Ue2=d2*U: positive down/equatorward, Eq(5.6) 
! un(3)=Ue3=d3*U: positive parallel to a field line, Eq(5.6) 
            DO jth=1,2
               Ue(i1d,jth) = Un_ms1(i,lp_plas,mp,jth)
            END DO

            DO jth=1,3
               
               IF (jth<=2) THEN
                  jth1=jth
                  jth2=jth
               ELSE IF (jth==3) THEN
                  jth1=1
                  jth2=2
               END IF
               dotprod = apexD(i,lp_plas,mp,east, jth1) * apexD(i,lp_plas,mp,east, jth2)  &
                    &  + apexD(i,lp_plas,mp,north,jth1) * apexD(i,lp_plas,mp,north,jth2)  &
                    &  + apexD(i,lp_plas,mp,up,   jth1) * apexD(i,lp_plas,mp,up   ,jth2)
               
               IF ( jth==1 ) THEN
                  apex_d1d1(i1d) = dotprod
               ELSE IF ( jth==2 ) THEN 
                  apex_d2d2(i1d) = dotprod
               ELSE IF ( jth==3 ) THEN 
                  apex_d1d2(i1d) = dotprod
               END IF
            END DO  !jth=1,2


! cross product of d1 & d2
! apex_D = ABS( D1 x D2 )  !eq(3.15)
!from apxntrpb4lf.f: line 1550
! E3(1) = D1(2)*D2(3) - D1(3)*D2(2)
!  d1xd2(1) = apexD(1,i,mp)%north * apexD(2,i,mp)%up    -  apexD(1,i,mp)%up    * apexD(2,i,mp)%north
! E3(2) = D1(3)*D2(1) - D1(1)*D2(3)
!  d1xd2(2) = apexD(1,i,mp)%up    * apexD(2,i,mp)%east  -  apexD(1,i,mp)%east  * apexD(2,i,mp)%up
! E3(3) = D1(1)*D2(2) - D1(2)*D2(1)
!  d1xd2(3) = apexD(1,i,mp)%east  * apexD(2,i,mp)%north -  apexD(1,i,mp)%north * apexD(2,i,mp)%east
! D = BHAT(1)*E3(1) + BHAT(2)*E3(2) + BHAT(3)*E3(3)
!???how can I get BHAT???
! BHAT(I) = B(I)/BMAG: unit vector along geomagnetic field direction
!  apex_D(i1d) = ABS( d1xd2(1)*d1xd2(1) + d1xd2(2)*d1xd2(2) + d1xd2(3)*d1xd2(3) )
!or, simply BE3 = BMAG/D
!  apex_D(i1d) = plasma_grid_3d(i,mp)%BM / Be3(1,mp,lp) 
            apex_D(i1d) = apexDscalar(i,lp_plas,mp)
            apex_BMAG(i1d) = plasma_grid_3d(i,lp_plas,mp,IBM)

            ni_oplus_1d(i1d)=plasma_3d(i,lp_plas,mp,1)
            ni_hplus_1d(i1d)=plasma_3d(i,lp_plas,mp,2)
            
            n_plus( i1d)=plasma_3d(i,lp_plas,mp,4)
            no_plus(i1d)=plasma_3d(i,lp_plas,mp,5)
            o2_plus(i1d)=plasma_3d(i,lp_plas,mp,6)
            n2_plus(i1d)=plasma_3d(i,lp_plas,mp,7)
            ti_oplus_1d(i1d)=plasma_3d(i,lp_plas,mp,11)
            
            tn(i1d)=TN_k(  i,lp_plas,mp) ![K]
            o(i1d)=ON_m3( i,lp_plas,mp) !m-3
            o2(i1d)=O2N_m3(i,lp_plas,mp)
            n2(i1d)=N2N_m3(i,lp_plas,mp)
            


         END DO i_loop

!nm20130830: Be3 is constant along magnetic field lines!
         Apex_BE3      =Be3(lp_plas,mp) 
!nm20130830:  Apex_BE3(CTIPDIM)=Be3(2,mp,lp_plas) !SH


         in1d = JMIN_IN(lp_plas) - JMIN_IN(lp_plas) + 1
         is1d = JMAX_IS(lp_plas) - JMIN_IN(lp_plas) + 1
         midpoint1d = ( JMAX_IS(lp_plas) - JMIN_IN(lp_plas) )/2 + 1
         NPTS = MaxFluxTube
         
         CALL calculate_field_line_integrals (NPTS,in1d,is1d,ds,midpoint1d, &
!nm20111107                                         Apex_d1,Apex_d2, &
     & Apex_BE3, &
                                         apex_D,apex_d1d1,apex_d1d2,apex_d2d2,apex_BMAG, &
                                         ni_oplus_1d,ni_hplus_1d,ti_oplus_1d,no_plus,o2_plus, &
                                         n2_plus,n_plus, &
                                         tn,o,o2,n2, &
!nm20111107                                         U_zonal,U_merid,U_vert, &
     & Ue, &
     & electron_density_out, &
                                         sigma_ped_1d,sigma_hall_1d, & !output
                                         sigma_phph_dsi_1d,sigma_lmlm_msi_1d, &
                                         sigma_h_1d,sigma_c_1d, &
                                         Kdmph_dsi_1d,Kdmlm_1d &
     &, Je_1d,mp,lp_plas )


! save to 3D array                                         
         call save2fli_array (lp_plas,mp, &
&                                         sigma_phph_dsi_1d,sigma_lmlm_msi_1d, &
&                                         sigma_h_1d,sigma_c_1d, &
&                                         Kdmph_dsi_1d,Kdmlm_1d )

!d!SMS$ignore begin
!d      print*,mype,'after save2fli_array'
!d!SMS$ignore end

! Tzu-Wei:save conductivities and Ue
          do i_set=in1d,is1d
                Ue1_3d(i_set,lp_dyn,mp) = Ue(i_set,1)
                Ue2_3d(i_set,lp_dyn,mp) = Ue(i_set,2)
                sigma_ped_3d(i_set,lp_dyn,mp) = sigma_ped_1d(i_set)
                sigma_hall_3d(i_set,lp_dyn,mp) = sigma_hall_1d(i_set)
                Ne_3d(i_set,lp_dyn,mp) = electron_density_out(i_set)
          enddo
!d!SMS$ignore begin
!d      print*,mype,lp_dyn,'after save 3d cond'
!d!SMS$ignore end

            EXIT lp_dyn_loop
!dbg20151107         END IF !imlat_dyn > imlat_plas ) THEN
      END DO      lp_dyn_loop!: DO lp_dyn=1,nmlat1

!d!SMS$ignore begin
!d      print*,mype,mp,'end of interface_field_line_integrals'
!d!SMS$ignore end

      END SUBROUTINE interface_field_line_integrals

      END MODULE module_interface_field_line_integrals
