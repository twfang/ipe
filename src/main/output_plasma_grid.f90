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
SUBROUTINE output_plasma_grid ( )
      USE module_precision
        USE module_IO,ONLY: filename
        USE module_FIELD_LINE_GRID_MKS,ONLY:plasma_grid_3d,plasma_grid_Z,plasma_grid_GL,JMIN_IN,JMAX_IS,mlon_rad,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,JMIN_IN_all,JMAX_IS_all
        USE module_IPE_dimension,ONLY: NPTS2D,NMP,NLP
IMPLICIT NONE
!-------------local
        CHARACTER (LEN=11) :: FORM_dum
        CHARACTER (LEN=7)  :: STATUS_dum
        INTEGER (KIND=int_prec) :: lun ,mp,lp,stat_alloc
      REAL (KIND=real_prec),DIMENSION(:,:), ALLOCATABLE :: dumm  !(NPTS2D,NMP)
!
      LUN=1006
      filename ='plasma_grid'
      FORM_dum ='unformatted' 
      STATUS_dum ='unknown'
!SMS$SERIAL BEGIN
      CALL open_file ( filename, LUN, FORM_dum, STATUS_dum ) 

    WRITE(UNIT=lun) NMP
    WRITE(UNIT=lun) NLP
    WRITE(UNIT=lun) NPTS2D
    WRITE(UNIT=lun) JMIN_IN(1:NLP)
    WRITE(UNIT=lun) JMAX_IS(1:NLP)
    WRITE(UNIT=lun) mlon_rad( 1: NMP+1 ) !rad
    WRITE(UNIT=lun) (plasma_grid_Z (JMIN_IN(lp):JMAX_IS(lp),lp),lp=1,NLP)  !meter
    WRITE(UNIT=lun) (plasma_grid_GL(JMIN_IN(lp):JMAX_IS(lp),lp),lp=1,NLP)  !rad
!SMS$SERIAL END


IF (.NOT.ALLOCATED(dumm) ) THEN
  ALLOCATE ( dumm(NPTS2D,NMP) &
     &,STAT=stat_alloc )         
      IF ( stat_alloc/=0 ) THEN
        print *,"sub-output_p:!STOP! ALLOCATION FAILD!:",stat_alloc
        STOP
      END IF
ELSE
STOP 'sub-output_p:!STOP! dumm has been allocated already???!!!'
END IF

!SMS$SERIAL(<plasma_grid_3d,IN> : default=ignore) BEGIN
do lp=1,NLP
  dumm(JMIN_IN_all(1,lp):JMAX_IS_all(1,LP),1:NMP) = plasma_grid_3d(JMIN_IN(lp):JMAX_IS(LP), lp,1:NMP,IGCOLAT)
enddo
WRITE(UNIT=lun) dumm !GCOLAT  !rad
do lp=1,NLP
  dumm(JMIN_IN_all(1,lp):JMAX_IS_all(1,LP),1:NMP) = plasma_grid_3d(JMIN_IN(lp):JMAX_IS(LP), lp,1:NMP,IGLON)
enddo
WRITE(UNIT=lun) dumm !GLON !rad
!SMS$SERIAL END

IF ( ALLOCATED(dumm) ) THEN 
  DEALLOCATE ( dumm &
     &,STAT=stat_alloc )         
      IF ( stat_alloc/=0 ) THEN
        print *,"sub-output_p:!STOP! DEALLOCATION FAILD!:",stat_alloc
        STOP
      END IF
ELSE
STOP 'sub-output_p:!STOP! dumm has not been allocated???!!!'
END IF


!SMS$SERIAL BEGIN
    CLOSE(UNIT=lun)
!SMS$SERIAL END
print *,'output_plasma_grid finished successfully!!!'
STOP
END SUBROUTINE output_plasma_grid
!---------------------------
