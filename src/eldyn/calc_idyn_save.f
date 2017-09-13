!This module calculates idyn_save which is the correspondence between lp & lp_dyn
!idyn_save is calculated on all processors but is only needed by the root in a serial region in module_plas2dyn_fli_array.f
      MODULE module_calc_idyn_save
      IMPLICIT NONE
      CONTAINS

      SUBROUTINE calc_idyn_save 
      USE module_precision
      USE cons_module               ,ONLY: xlatm,idyn_save,lp_dyn_eq
      USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_grid_GL,JMAX_IS
      USE module_physical_constants ,ONLY: rtd
      USE module_input_parameters   ,ONLY: mpstop,mype
      USE module_IPE_dimension      ,ONLY: NMP,NLP

      REAL    (KIND=real_prec) :: mlat_dyn, mlat_plas
      INTEGER (KIND=int_prec ) :: mp,lp,lp_dyn,imlat_dyn,imlat_plas

      idyn_save  = 0
      DO mp = 1,mpstop
        DO lp = 1,NLP
          mlat_plas  = 90. - plasma_grid_GL(JMAX_IS(lp),lp)*rtd
          imlat_plas = INT(mlat_plas*10.)
          lp_dyn_loop: DO lp_dyn=1,lp_dyn_eq !from SH toward eq
            mlat_dyn  = xlatm(lp_dyn)*rtd  ![deg]
            imlat_dyn = INT(mlat_dyn*10.)
            IF ( lp < NLP ) THEN
              IF ( imlat_dyn > imlat_plas ) THEN
                EXIT lp_dyn_loop
              ELSEIF ( imlat_dyn < imlat_plas ) THEN
                CYCLE lp_dyn_loop
              ENDIF
              idyn_save(lp_dyn)=lp  !correspondance between lp & lp_dyn
            ELSEIF ( lp == nlp ) THEN
              IF ( lp_dyn < lp_dyn_eq ) THEN
                CYCLE lp_dyn_loop
              ENDIF
              idyn_save(lp_dyn)=lp  !correspondance between lp & lp_dyn
            ENDIF
          ENDDO lp_dyn_loop
        ENDDO
      ENDDO

      END SUBROUTINE calc_idyn_save

      END MODULE module_calc_idyn_save
