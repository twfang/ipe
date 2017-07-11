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
      module module_save2fli_array
        IMPLICIT NONE
        PRIVATE
        PUBLIC :: save2fli_array
      contains
        subroutine save2fli_array (lp_plas,mp, &
     &                                         sigma_phph_dsi_1d,sigma_lmlm_msi_1d, &
     &                                         sigma_h_1d,sigma_c_1d, &
     &                                         Kdmph_dsi_1d,Kdmlm_1d )
          USE module_precision
          USE module_eldyn,ONLY: plas_fli !t,Je_3d
          use module_input_parameters,ONLY:mype 
          IMPLICIT NONE
          INTEGER (KIND=int_prec), intent(in)  ::  mp
          INTEGER (KIND=int_prec), intent(in)  ::  lp_plas
!
          REAL (KIND=real_prec), intent(in)  ::    sigma_phph_dsi_1d(2)      !(5.13) divided by |sin I_m |
          REAL (KIND=real_prec), intent(in)  ::    sigma_lmlm_msi_1d(2)      !(5.14) multiplied by | sin I_m |
          REAL (KIND=real_prec), intent(in)  ::	   sigma_h_1d(2)      !(5.17)
          REAL (KIND=real_prec), intent(in)  ::	   sigma_c_1d(2)      !(5.18)
          REAL (KIND=real_prec), intent(in)  ::    Kdmph_dsi_1d(2)      !(5.19) divided by |sin I_m |
          REAL (KIND=real_prec), intent(in)  ::	   Kdmlm_1d(2)	  !(5.20) plus or minus ????
!
          INTEGER (KIND=int_prec) ::  ihem
!
!SH ihem=2; NH ihem=1
          do ihem=1,2
                plas_fli(ihem,lp_plas,mp,1) = sigma_phph_dsi_1d(ihem)
                plas_fli(ihem,lp_plas,mp,2) = sigma_lmlm_msi_1d(ihem)
                plas_fli(ihem,lp_plas,mp,3) = sigma_h_1d(       ihem)
                plas_fli(ihem,lp_plas,mp,4) = sigma_c_1d(       ihem)
                plas_fli(ihem,lp_plas,mp,5) = Kdmph_dsi_1d(     ihem)
                plas_fli(ihem,lp_plas,mp,6) = Kdmlm_1d(         ihem)
          end do
! Tzu-Wei TEST
!SMS$IGNORE begin
      print *,mype,ihem,lp_plas,mp,'TEST zigm within plasma folder'
      print *,mype,'zigm11',MAXVAL(plas_fli(:,:,:,1)),MINVAL(plas_fli(:,:,:,1))
      print *,mype,'zigm22',MAXVAL(plas_fli(:,:,:,2)),MINVAL(plas_fli(:,:,:,2))
      print *,mype,'zigm2',MAXVAL(plas_fli(:,:,:,3)),MINVAL(plas_fli(:,:,:,3))
      print *,mype,'zigmc',MAXVAL(plas_fli(:,:,:,4)),MINVAL(plas_fli(:,:,:,4))
      print *,mype,'rim1',MAXVAL(plas_fli(:,:,:,5)),MINVAL(plas_fli(:,:,:,5))
      print *,mype,'rim2',MAXVAL(plas_fli(:,:,:,6)),MINVAL(plas_fli(:,:,:,6))
!SMS$IGNORE end

!t          IF ( sw_3DJ==1 ) THEN
!t            DO jth=1,2
!t              Je_3d(IN:IS,mp,lp_dyn,jth) = Je_1d(1:CTIPDIM,jth)
!t            END DO
!t          END IF !( sw_3DJ==1 ) THEN 
!t          DO i=JMIN_IN(lp_plas),JMAX_IS(lp_plas)
!t            i1d=i-JMIN_IN(lp_plas)+1
!t            sigma_ped(i,lp_dyn,mp) =sigma_ped_1d(i1d)                              
!t            sigma_hall(i,lp_dyn,mp)=sigma_hall_1d(i1d)                              
!t          END DO

        end subroutine save2fli_array
      end module module_save2fli_array
