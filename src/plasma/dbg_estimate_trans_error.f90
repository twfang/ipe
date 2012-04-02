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
!NOTE: dbg20120228: only temporary used to estimate the error created by the perp transport for debug purpose
      subroutine dbg_estimate_trans_error ( utime )
        USE module_precision
        USE module_IPE_dimension,ONLY:NMP_all,NLP_all,IPDIM
        USE module_PLASMA,ONLY:plasma_3d
        IMPLICIT NONE
        INTEGER (KIND=int_prec),  INTENT(IN)  :: utime    !universal time [sec]

!---local variables---
        INTEGER (KIND=int_prec) :: lp,mp
        REAL (KIND=real_prec),dimension(NMP_all,NLP_all) :: opmaxpc,opminpc
        REAL (KIND=real_prec),PARAMETER :: fixedval=100.0

!
do lp=1,NLP_all
do mp=1,NMP_all
opmaxpc(mp,lp)=( &
     MAXVAL( plasma_3d(mp,lp)%N_m3( 1,1:IPDIM) )-fixedval &
              & )/fixedval
opminpc(mp,lp)=( &
     MINVAL( plasma_3d(mp,lp)%N_m3( 1,1:IPDIM) )-fixedval &
              & )/fixedval
enddo
enddo

write(5000,fmt=*) utime,MAXVAL(opmaxpc),MINVAL(opminpc)


END       subroutine dbg_estimate_trans_error 
