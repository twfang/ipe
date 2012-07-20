!20110911: note: openmp was tried on jet but did not work: only thread 0 was used not the other thread...although other threads did exist...needs more investigation...
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
      MODULE module_PLASMA
      USE module_precision
      USE module_IPE_dimension,ONLY: ISPEC,ISPET,ISPEV,IPDIM, NPTS2D,NLP
      IMPLICIT NONE
! --- PRIVATE ---
!
! --- PUBLIC ---
!solver REALL needs these parameters from the previous time step
      !.. TE_TI_k(3,J) = Te, TE_TI_k(2,J) = Ti = TE_TI_k(2,J) [kelvin]
      TYPE :: plasma_data_1d
        REAL(KIND=real_prec), DIMENSION(      IPDIM) :: Te_k  !TE_TI(3)
        REAL(KIND=real_prec), DIMENSION(ISPET,IPDIM) :: Ti_k  !TE_TI(1:2)
        REAL(KIND=real_prec), DIMENSION(ISPEC,IPDIM) :: N_m3
!dbg20110927        REAL(KIND=real_prec), DIMENSION(ISPEV,IPDIM) :: V_ms1
!dbg20110927        REAL(KIND=real_prec), DIMENSION(      IPDIM) :: heating_rate_e_cgs ![eV cm-3 s-1] !EHT(3)
!dbg20110927        REAL(KIND=real_prec), DIMENSION(ISPET,IPDIM) :: heating_rate_i_cgs  !EHT(1:2)
!dbg20110923        REAL(KIND=real_prec), DIMENSION(      IPDIM) :: NO_m3 
      END TYPE  plasma_data_1d
      TYPE(plasma_data_1d),ALLOCATABLE,TARGET,PUBLIC :: plasma_3d(:,:,:)

!neutral needs these parameters from the ionosphere in addition to N&T
      TYPE :: plasma_data_1d4n
        REAL(KIND=real_prec), DIMENSION(ISPEV) :: V_ms1
!???      REAL(KIND=real_prec) :: NHEAT
      END TYPE  plasma_data_1d4n
      TYPE(plasma_data_1d4n),ALLOCATABLE,TARGET,PUBLIC :: plasma_3d4n(:,:) !(NPTS2D, NMP)

      TYPE(plasma_data_1d), PUBLIC :: n0_1d !N&T after perpendicular transport

!only for debug, o+
!d      REAL(KIND=real_prec), DIMENSION(NPTS2D), PUBLIC :: n0_2dbg 
!      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP), PUBLIC :: TE_TI_k
!      REAL(KIND=real_prec), DIMENSION(ISPEC,NPTS2D,NMP), PUBLIC :: XIONN_m3
!      REAL(KIND=real_prec), DIMENSION(ISPEC,NPTS2D,NMP), PUBLIC :: XIONV_ms1

!V_ExB m/s at the apex height
      REAL(KIND=real_prec), DIMENSION(:,:),ALLOCATABLE, PUBLIC :: VEXBup !DIMENSION(NMP,NLP)

! save ut so that other subroutines can refer to it
      INTEGER (KIND=int_prec),PUBLIC:: utime_save

      PRIVATE
      PUBLIC :: plasma,plasma_data_1d,plasma_data_1d4n


      CONTAINS
!---------------------------
      SUBROUTINE plasma ( utime )
      USE module_input_parameters,ONLY:lpstrt,lpstop,lpstep,mpstrt,mpstop,mpstep,ip_freq_output,start_time,stop_time,sw_neutral_heating_flip,sw_perp_transport,lpmin_perp_trans,lpmax_perp_trans,sw_para_transport,sw_debug &
&, sw_dbg_perp_trans , sw_exb_up
      USE module_heating_rate,ONLY: hrate_cgs_save
      USE module_physical_constants,ONLY:rtd,zero
      USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,plasma_grid_3d,plasma_grid_GL,mp_save,lp_save,plasma_grid_Z,JMAX_IS
      IMPLICIT NONE
!------------------------
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
!--- local variables ---
      INTEGER (KIND=int_prec) :: mp
      INTEGER (KIND=int_prec) :: lp
      INTEGER (KIND=int_prec) :: i,j,midpoint
!d      INTEGER :: lun_dbg=999
!t      REAL(KIND=real_prec) :: phi_t0   !magnetic longitude,phi at T0
!t      REAL(KIND=real_prec) :: theta_t0 !magnetic latitude,theta at T0
!

!---------------
! array initialization
!dbg20120313 note! needed to comment out temporarily for exb reading test:sw_exb=5
if ( sw_exb_up/=5 ) then
      VEXBup(:,:)=zero
end if


! save ut so that other subroutines can refer to it
      utime_save=utime

      apex_longitude_loop: DO mp = mpstrt,mpstop,mpstep !1,NMP
        mp_save=mp
        IF ( sw_neutral_heating_flip==1 ) hrate_cgs_save=zero
        if ( sw_debug )  WRITE (0,"('sub-p: mp=',I4)")mp
!d        n0_2dbg(:)=zero

        DO i=1,NPTS2D
          plasma_3d4n(i,mp)%V_ms1(1:ISPEV)=zero
        END DO


!!!dbg20120125: only temporary used to switch on the transport only during the daytime...
!JFM    IF ( sw_rw_sw_perp_trans.AND.sw_perp_transport(mp)==0 )  CALL activate_perp_transport (utime,mp)
!!!dbg20120125:

        apex_latitude_height_loop: DO lp = lpstrt,lpstop,lpstep
          lp_save=lp


!dbg20120228: debug how2validate the transport
if(sw_dbg_perp_trans.and.utime==start_time.and.lp==lpstrt)then
DO j=1,NLP
DO i=1,IPDIM
plasma_3d(mp,j)%N_m3( 1:ISPEC,i)=100.0
plasma_3d(mp,j)%Te_k(         i)=100.0
plasma_3d(mp,j)%Ti_k( 1:ISPET,i)=100.0
END DO !i=1,IPDIM
END DO !j=1,NLP
end if
!if(sw_dbg_perp_trans) print *, '1!dbg max o+',MAXVAL( plasma_3d(mp,lp)%N_m3( 1,1:IPDIM) ),MINVAL( plasma_3d(mp,lp)%N_m3( 1,1:IPDIM) )




!20111025: not sure if these lines work when ut=0 & HPEQ=0.5(initial profiles are prepared within flip) , or maybe it is ok if they are zero?
!save the values from the previous time step...
          DO i=1,IPDIM
            n0_1d%N_m3( 1:ISPEC,i)=plasma_3d(mp,lp)%N_m3( 1:ISPEC,i)
            n0_1d%Te_k(         i)=plasma_3d(mp,lp)%Te_k(         i)
            n0_1d%Ti_k( 1:ISPET,i)=plasma_3d(mp,lp)%Ti_k( 1:ISPET,i)
          END DO

!if(sw_dbg_perp_trans) print *, '2!dbg max o+',MAXVAL( n0_1d%N_m3( 1,1:IPDIM) ),MINVAL( n0_1d%N_m3( 1,1:IPDIM) )

          if ( sw_debug )  WRITE (0,"('sub-p: lp=',I4)")lp





          IF ( sw_perp_transport(mp)>=1 ) THEN
            IF ( lp>=lpmin_perp_trans.AND.lp<=lpmax_perp_trans ) THEN
              CALL perpendicular_transport ( utime,mp,lp )

!dbg20110927:put the values back to the 3D array!
!dbg plasma_3d(mp,lp)%N_m3(1:ISPEC,1:IPDIM)=n0_1d%N_m3(1:ISPEC,1:IPDIM)


            ELSE  !IF ( lp>lpmin_perp_trans ) THEN
if(utime==start_time) then
midpoint=JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
print "('NO PERP. TRANS: mp=',I3,' lp=',I4,' mlatNd',F8.3,' apht',F8.2)", mp,lp,(90.-plasma_grid_GL(JMIN_IN(lp),lp)*rtd) &
& ,plasma_grid_Z(midpoint,lp)*1.0e-3
endif
            END IF
          END IF !( sw_perp_transport>=1 ) THEN
          

! update the boundary conditions if the top of the flux tube is open
!t        CALL update_flux_tube_boundary_condition ( )

!!nm20110822: moved to flux_tube_plasma  
!! calculate Solar Zenith Angle [radians]
!!          CALL Get_SZA ( utime,mp,lp ) 

! call flux tube solver
          IF ( sw_para_transport==1 ) THEN 

            CALL flux_tube_solver ( utime,mp,lp )
          ELSE IF ( sw_para_transport==0 ) THEN 

!dbg20111101:v9: temporary ...
            DO i=1,IPDIM
              plasma_3d(mp,lp)%N_m3( 1:ISPEC,i)=            n0_1d%N_m3( 1:ISPEC,i)
              plasma_3d(mp,lp)%Te_k(         i)=            n0_1d%Te_k(         i)
              plasma_3d(mp,lp)%Ti_k( 1:ISPET,i)=            n0_1d%Ti_k( 1:ISPET,i)
            END DO
          END IF !( sw_para_transport==1 ) THEN           

! calculate neutral heating rate: NHEAT_mks in [eV kg-1 s-1]
!20110729: temporarily commented out to save time...
!dbg20110927          IF ( sw_neutral_heating_flip==1 ) &
!dbg20110927     &       CALL get_neutral_heating_rate ( )

! calculate the field line integrals for the electrodynamic solver
!t        CALL calculate_field_line_integrals ( )

        END DO apex_latitude_height_loop !: DO lp = 1
      END DO apex_longitude_loop !: DO mp = 

!dbg20120228: debug how2validate the transport
if(sw_dbg_perp_trans) call dbg_estimate_trans_error (utime)

! output plasma parameters to a file
      IF ( MOD( (utime-start_time),ip_freq_output)==0 ) THEN 
if(sw_debug) print *,'before call to output plasma',utime,start_time,ip_freq_output
!dbg20110923segmentation fault??? memory allocation run time error???
        CALL io_plasma_bin ( 1, utime )

!dbg20110927: o+ only
!d IF ( sw_perp_transport>=1 ) THEN
!d if (utime==start_time)  open(unit=lun_dbg,file='dbg_trans',status='unknown',form='formatted')
!d write(unit=lun_dbg,fmt='(20E12.4)') n0_2dbg(1:npts2d)
!d if (utime==stop_time)   close(unit=lun_dbg)
!d END IF !( sw_perp_transport>=1 ) THEN

      END IF      !IF ( MOD( (utime-start_time),ip_freq_output)==0 ) THEN 

      END SUBROUTINE plasma
      END MODULE module_PLASMA
