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
      SUBROUTINE flux_tube_solver ( utime,mp,lp )
      USE module_precision
      USE module_IPE_dimension,ONLY: ISPEC,ISPEV,IPDIM
      USE module_FIELD_LINE_GRID_MKS,ONLY: JMIN_IN,JMAX_IS,plasma_grid_3d,plasma_grid_Z,plasma_grid_GL,Pvalue
      USE module_input_parameters,ONLY: time_step,F107D,F107AV,DTMIN_flip  &
     &, sw_INNO,FPAS_flip,HPEQ_flip,HEPRAT_flip,COLFAC_flip,sw_IHEPLS,sw_INPLS,sw_debug,iout, start_time, sw_wind_flip, sw_depleted_flip, start_time_depleted, sw_output_fort167
      USE module_NEUTRAL_MKS,ONLY: ON_m3,HN_m3,N2N_m3,O2N_m3,HE_m3,N4S_m3 &
     &, TN_k,TINF_k,un_ms1
      USE module_PLASMA,ONLY: plasma_3d, plasma_1d !dbg20120501
!dbg20110927      USE module_heating_rate,ONLY: NHEAT_cgs
      USE module_physical_constants,ONLY: pi,zero
      USE module_IO,ONLY: PRUNIT,LUN_FLIP1,LUN_FLIP2,LUN_FLIP3,LUN_FLIP4
      USE module_unit_conversion,ONLY: M_TO_KM

      IMPLICIT NONE
!------------------------
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      INTEGER (KIND=int_prec), INTENT(IN) :: mp    !longitude
      INTEGER (KIND=int_prec), INTENT(IN) :: lp    !latitude
      REAL (KIND=real_prec) :: ltime !local time [hour]

!--- for CTIPINT
      INTEGER, POINTER :: IN,IS !.. lcv + spatial grid indices
      INTEGER JMINX,JMAXX !.. lcv + spatial grid indices
      INTEGER CTIPDIM         !.. CTIPe array dimension, must equal to FLDIM
      DOUBLE PRECISION ::  PCO
      INTEGER ::  INNO            !.. switch to turn on FLIP NO calculation if <0
      DOUBLE PRECISION, DIMENSION(IPDIM) ::  ZX,GLX,SLX,BMX,GRX,OX,HX,N2X,O2X,HEX,N4SX,NNOX &
     &, UNX &     ! assumed a component parallel to a field line [m s-1]
     &, TNX &
       !.. TINFX has to be an array for grazing incidence column densities
     &, TINFX & !.. Exospheric temperature [k]
     &, SZA_dum
      DOUBLE PRECISION :: DTMIN !.. Minimum time step allowed (>=10 secs?)
      DOUBLE PRECISION DT

      REAL ::  F107D_dum,F107A_dum         !.. daily and 81 day average F10.7      !.. 
      DOUBLE PRECISION ::  FPAS,HEPRAT,COLFACX
      DOUBLE PRECISION HPEQ

      INTEGER ::  IHEPLS,INPLS  !.. switches He+ and N+ diffusive solutions on 
      DOUBLE PRECISION &
      !.. EHTX(3,J) = e heating rate, EHTX(1,J) = ion heating rate, EHTX(2,J) unused
     &  EHTX(3,IPDIM) &
      !.. TE_TI(3,J) = Te, TE_TIX(2,J) = Ti = TE_TIX(2,J)
     & ,TE_TIX(3,IPDIM) &
     & ,XIONNX(ISPEC,IPDIM),XIONVX(ISPEC,IPDIM) &
     & ,NHEAT(IPDIM)   !.. Neutral heating rate [eV/cm^3/s]

      INTEGER EFLAG(11,11)    !.. error flags, check =0 on return from FLIP
      INTEGER :: PRUNIT_dum !.. Unit number to print results
      INTEGER (KIND=int_prec) :: midpoint
!      INTEGER (KIND=int_prec) :: stat_alloc
      INTEGER (KIND=int_prec) :: ipts,i
      INTEGER (KIND=int_prec),PARAMETER :: ip_freq_output_fort=900      
      INTEGER (KIND=int_prec) :: jth !dbg20120501
!----------------------------------





! frequency of output to fort167/8 170/1 within CTIP-int
!d      IF ( lp==170.and.mp==1.and.MOD( (utime-start_time), ip_freq_output_fort)==0 ) THEN
!d        sw_output_fort167=.TRUE. 
!d      ELSE 
!d        sw_output_fort167=.FALSE.
!d      END IF
 
      IN => JMIN_IN(lp)
      IS => JMAX_IS(lp)

! make sure that JMINX equals to 1
      JMINX   = IN - IN + 1
      JMAXX   = IS - IN + 1
      CTIPDIM = IS - IN + 1   

!nm20110822: no more allocatable arrays
!      ALLOCATE ( ZX(JMINX:JMAXX),  SLX(JMINX:JMAXX),  GLX(JMINX:JMAXX), BMX(JMINX:JMAXX), GRX(JMINX:JMAXX) &
!     &          ,OX(JMINX:JMAXX),   HX(JMINX:JMAXX),  N2X(JMINX:JMAXX), O2X(JMINX:JMAXX), HEX(JMINX:JMAXX) &
!     &          ,N4SX(JMINX:JMAXX),TNX(JMINX:JMAXX),TINFX(JMINX:JMAXX), UNX(JMINX:JMAXX), SZA_dum(JMINX:JMAXX) &
!     &          ,NNOX(JMINX:JMAXX),EHTX(1:3,JMINX:JMAXX),TE_TIX(1:3,JMINX:JMAXX) &
!     &          ,XIONNX(1:ISPEC,JMINX:JMAXX),XIONVX(1:ISPEC,JMINX:JMAXX),NHEAT(JMINX:JMAXX)  &
!     & ,STAT=stat_alloc         )
!if ( stat_alloc/=0 ) then
!  print *, ALLOCATED( ZX )
!  print *,"!STOP! ALLOCATION FAILD! in flux_tube_solver:",stat_alloc,mp,lp,in,is,jminx,jmaxx,ctipdim
!  STOP
!endif


      
      ZX(1:CTIPDIM)  = plasma_grid_Z(IN:IS) * M_TO_KM !convert from m to km 
      PCO = Pvalue(lp)  !Pvalue is a single value
      SLX(1:CTIPDIM) = plasma_grid_3d(IN:IS,mp)%SL
      GLX(1:CTIPDIM) = pi/2. - plasma_grid_GL(IN:IS)  ! magnetic latitude [radians]
      BMX(1:CTIPDIM) = plasma_grid_3d(IN:IS,mp)%BM   !Tesla
      GRX(1:CTIPDIM) = plasma_grid_3d(IN:IS,mp)%GR

      OX(1:CTIPDIM) = ON_m3(IN:IS,mp) !(m-3)
      HX(1:CTIPDIM) = HN_m3(IN:IS,mp)
      N2X(1:CTIPDIM) = N2N_m3(IN:IS,mp)
      O2X(1:CTIPDIM) = O2N_m3(IN:IS,mp)
      HEX(1:CTIPDIM) = HE_m3(IN:IS,mp)
      N4SX(1:CTIPDIM) = N4S_m3(IN:IS,mp)
      
      INNO = sw_INNO

      TNX(1:CTIPDIM) = TN_k(IN:IS,mp)
      TINFX(1:CTIPDIM) = TINF_k(IN:IS,mp)


! FLIP assumes positive SOUTHWARD along a field line
      IF ( sw_wind_flip == 1 ) THEN
        UNX(1:CTIPDIM)  = (-1.) * Un_ms1(3,IN:IS,mp) 
      ELSE IF ( sw_wind_flip == 0 ) THEN
        UNX(1:CTIPDIM)  = 0.0
      END IF

      DT        = REAL(time_step)
      DTMIN     = DTMIN_flip
      F107D_dum = F107D
      F107A_dum = F107AV

!nm20110822: moved from module_plasma
!! I need to get the new get_sza based on phil's method
!! calculate Solar Zenith Angle [radians]
      CALL Get_SZA ( utime,mp,lp, SZA_dum )
!      SZA_dum(JMINX:JMAXX)   = SZA_rad(IN:IS)
!nm20110822: no more allocatable arrays
!      IF ( ALLOCATED( SZA_rad  ) )  DEALLOCATE( SZA_rad, STAT=stat_alloc  )
!if ( stat_alloc/=0 ) then
!  print *, ALLOCATED( sza_rad )
!  print *,"!STOP! SZA_rad DEALLOCATION FAILED! in flux_tube_solver:",stat_alloc,mp,lp,in,is,jminx,jmaxx
!  STOP
!endif

      FPAS      = FPAS_flip

!nm110510: test the depleted flux tube
      IF ( utime == start_time ) THEN
        HPEQ      = HPEQ_flip
      ELSE
        HPEQ      = 0.0

        IF ( sw_depleted_flip==1 .AND. utime == start_time_depleted ) THEN
          HPEQ      = - 0.1
        ENDIF
      ENDIF
      

      HEPRAT    = HEPRAT_flip
      COLFACX   = COLFAC_flip

      IHEPLS    = sw_IHEPLS
      INPLS     = sw_INPLS
  
!! array initialization
      EFLAG(:,:)=0

!dbg20110802: 3D multiple-lp run
!if (mp==1 .and.lp==10) then
!  sw_debug=.true.
!else
!  sw_debug=.false.
!end if


IF ( sw_debug ) THEN

print *,'sub-flux_tube_solver'
print "('mp=',i6,' lp=',i6,' JMINX=',I6,' JMAXX=',I6)", mp,lp,jminx,jmaxx
print "('CTIPDIM=',I6)", ctipdim
print "('Z [km]     =',2F10.4)", ZX(jminx),ZX(jmaxx)
print "('PCO        =',2F10.4)", PCO,Pvalue(lp)
print "('SLX [m]    =',2E12.4)", SLX(jminx), SLX(jmaxx) 
print "('GLX [deg]  =',2F10.4)",(GLX(jminx)*180./pi),(GLX(jmaxx)*180./pi)
print "('BMX [Tesla]    =',2E12.4)", BMX(jminx), BMX(jmaxx)
print "('GRX[m2 s-1]=',2E12.4)",GRX(jminx),GRX(jmaxx)
!---neutral parameters
print "('LOG10 OX [m-3]     =',2F10.4)",LOG10(OX(jminx)),LOG10(OX(jmaxx))
print "('LOG10 HX [m-3]     =',2F10.4)",LOG10(HX(jminx)),LOG10(HX(jmaxx))
print "('LOG10 N2X [m-3]    =',2F10.4)",LOG10(N2X(jminx)),LOG10(N2X(jmaxx))
print "('LOG10 O2X [m-3]    =',2F10.4)",LOG10(O2X(jminx)),LOG10(O2X(jmaxx))
print "('LOG10 HEX [m-3]    =',2F10.4)",LOG10(HEX(jminx)),LOG10(HEX(jmaxx))
print "('LOG10 N4SX [m-3]   =',2F10.4)",LOG10(N4SX(jminx)),LOG10(N4SX(jmaxx))

print "('INNO =',I6)",INNO
IF ( INNO>=0 )  & !when CTIPe calculates NO
     & print "('LOG10 NNOX [m-3]   =',2F10.4)",LOG10(NNOX(jminx)),LOG10(NNOX(jmaxx))
print "('Tn [K]       =',2F10.4)",TNX(jminx),TNX(jmaxx)
print "('TINF [K]     =',2F10.4)",TINFX(jminx),TINFX(jmaxx)
print "('UNX [m s-1]  =',2F10.4)",UNX(jminx),UNX(jmaxx)

print "('DT [sec]     =',F10.4)",DT
print "('DTMIN [sec]  =',F10.4)",DTMIN
print "('F107D_dum    =',F10.4)",F107D_dum
print "('F107A_dum    =',F10.4)",F107A_dum
print "('SZA [deg]    =',2F10.4)",SZA_dum(jminx)*180./pi,SZA_dum(jmaxx)*180./pi
print "('FPAS         =',F10.4)",FPAS
print "('HPEQ         =',F10.4)",HPEQ
print "('HEPRAT       =',F10.4)",HEPRAT
print "('COLFACX      =',F10.4)",COLFACX
print "('IHEPLS       =',I6)",IHEPLS
print "('INPLS        =',I6)",INPLS

END IF !( sw_debug ) then


IF( sw_output_fort167 ) then
!dbg20120125:      midpoint = JMINX + (CTIPDIM-1)/2
      midpoint = IN + (IS-IN)/2
      IF ( lp>=1 .AND. lp<=6 )  midpoint = midpoint - 1
!nm20110909: calculating LT should be made FUNCTION!!!
      ltime = REAL(utime)/3600.0 + (plasma_grid_3d(midpoint,mp)%GLON*180.0/pi)/15.0
      IF ( ltime > 24.0 )  ltime = MOD(ltime, 24.0)

      WRITE(UNIT=LUN_FLIP1,FMT="('mp=',i3,' lp=',i3,' U',i3,' North, UT=',2F10.3)") mp,lp,LUN_FLIP1,REAL(UTIME)/3600., ltime
      WRITE(UNIT=LUN_FLIP2,FMT="('mp=',i3,' lp=',i3,' U',i3,' North, UT=',2F10.3)") mp,lp,LUN_FLIP2,REAL(UTIME)/3600., ltime
      WRITE(UNIT=LUN_FLIP3,FMT="('mp=',i3,' lp=',i3,' U',i3,' South, UT=',2F10.3)") mp,lp,LUN_FLIP3,REAL(UTIME)/3600., ltime
      WRITE(UNIT=LUN_FLIP4,FMT="('mp=',i3,' lp=',i3,' U',i3,' South, UT=',2F10.3)") mp,lp,LUN_FLIP4,REAL(UTIME)/3600., ltime

print *,'sub-fl: UTs=',UTIME,' LThr=',ltime,' mp',mp,' lp',lp
END IF !( sw_debug ) then



!dbg20110131:
IF ( sw_debug )  WRITE(UNIT=PRUNIT,FMT="('mp=',i6,' lp=',i6,' UT=',F10.2)") mp,lp,REAL(UTIME)/3600.

      DO ipts=1,CTIPDIM
!N&T from the previous time step are absolute necesary for the solver...
!dbg20120501
        DO jth=1,ISPEC
          XIONNX(jth,ipts) = plasma_1d(jth,ipts)
        END DO !jth
!te
        TE_TIX(3,ipts) = plasma_1d(ISPEC+1,ipts)
!ti
        DO jth=1,2
          TE_TIX(jth,ipts) = plasma_1d(jth+ISPEC+1,ipts)
        END DO !jth
!vi
        DO jth=1,ISPEC
IF ( jth<=2 ) THEN
          XIONVX(jth,ipts) = plasma_1d(jth+ISPEC+3,ipts) 
ELSE
          XIONVX(jth,ipts) = zero
END IF
        END DO !jth
!dbg20120501         XIONNX(1:ISPEC,ipts) = n0_1d%N_m3(1:ISPEC,ipts)
!dbg20120501         TE_TIX(3      ,ipts) = n0_1d%Te_k(        ipts)
!dbg20120501         TE_TIX(2      ,ipts) = n0_1d%Ti_k(      2,ipts)
!dbg20120501         TE_TIX(1      ,ipts) = n0_1d%Ti_k(      1,ipts) 

!### need to change these derived data type to a simpler arrays!!!
!nm20120412: need to restore the save V//(1:2) for parallel conv. effect
!nm20120412         XIONVX(1:ISPEC,ipts) = zero  !dbg20110927
!dbg20120501: add v// to perp transport
!dbg20120501         XIONVX(1,ipts)=plasma_1d(1,1,ipts) 
!dbg20120501         XIONVX(2,ipts)=plasma_1d(2,1,ipts) 
!dbg20120501         XIONVX(3:ISPEC,ipts) = zero

!dbg20110927         EHTX(3          ,ipts)= plasma_3d(mp,lp)%heating_rate_e_cgs(  ipts)
!dbg20110927         EHTX(2          ,ipts)= plasma_3d(mp,lp)%heating_rate_i_cgs(2,ipts)
!dbg20110927         EHTX(1          ,ipts)= plasma_3d(mp,lp)%heating_rate_i_cgs(1,ipts)

!20110930: inclusion for future neutral coupling:
! auroral electron heating-->EHTX(3)
! frictional heating for ions-->EHTX(1) 

         EHTX(  1:3        ,ipts)=zero  !dbg20110927

!dbg20110809 v2
!dbg20110923     NNOX(          ipts)= plasma_3d(mp,lp)%NO_m3(        ipts)
         IF ( INNO<0 ) THEN   !when flip calculates NO
           NNOX(              ipts)=zero !dbg20110927
         ELSE !when CTIPe calculates NO
           print *,'CTIPe calculates NO'
         END IF
         NHEAT(             ipts)=zero !dbg20110927
      END DO !ipts=

! call the flux tube solver (FLIP)
      CALL CTIPINT( &
     &             JMINX, & !.. index of the first point on the field line
     &             JMAXX, & !.. index of the last point on the field line
     &           CTIPDIM, & !.. CTIPe array dimension, must equal to FLDIM
     &                ZX(1:CTIPDIM), & !.. array, altitude (km)
     &               PCO, & !.. p coordinate (L-shell)
     &               SLX(1:CTIPDIM), & !.. array, distance of point from northern hemisphere (meter)
     &               GLX(1:CTIPDIM), & !.. array, magnetic latitude (radians)
     &               BMX(1:CTIPDIM), & !.. array, magnetic field strength, (Tesla)
     &               GRX(1:CTIPDIM), & !.. array, gravity, m2 s-1
     &                OX(1:CTIPDIM), & !.. array, O density (m-3)
     &                HX(1:CTIPDIM), & !.. array, H density (m-3)
     &               N2X(1:CTIPDIM), & !.. array, N2 density (cm-3)
     &               O2X(1:CTIPDIM), & !.. array, O2 density (cm-3)
     &               HEX(1:CTIPDIM), & !.. array, He density (cm-3)
     &              N4SX(1:CTIPDIM), & !.. array, N(4S) density (cm-3)
     &              INNO, & !.. switch to turn on FLIP NO calculation if <0
     &              NNOX(1:CTIPDIM), & !.. array, NO density (cm-3)
     &               TNX(1:CTIPDIM), & !.. array, Neutral temperature (K)
     &             TINFX(1:CTIPDIM), & !.. array, Exospheric Neutral temperature (K)
     &               UNX(1:CTIPDIM), & !.. array, Neutral wind (m/s), positive northward,horizontal component in the magnetic meridian??? or component parallel to a field line???
     &                DT, & !.. CTIPe time step (secs)
     &             DTMIN, & !.. Minimum time step allowed (>=10 secs?)
     &         F107D_dum, & !.. Daily F10.7
     &         F107A_dum, & !.. 81 day average F10.7
     &           SZA_dum(1:CTIPDIM), & !.. Solar Zenith angle (radians)
     &              FPAS, & !.. Pitch angle scattering fraction
     &              HPEQ, & !.. Sets initial equatorial H+ density. See declaration below
     &            HEPRAT, & !.. Intial He+/H+ ratio (.01 to 1.0)
     &           COLFACX, & !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)
     &            IHEPLS, & !.. switches He+ diffusive solution on if > 0
     &             INPLS, & !.. switches N+ diffusive solution on if > 0
     &              EHTX(1:3,1:CTIPDIM), & !.. IN/OUT 2D array, Electron & ion heating rate (eV cm-3 s-1)
     &            TE_TIX(1:3,1:CTIPDIM), & !.. IN/OUT: 2D array, Electron and ion temperatures (K) (see below)
     &     XIONNX(1:ISPEC,1:CTIPDIM),XIONVX(1:ISPEC,1:CTIPDIM), & !.. IN/OUT: 2D array, Storage for ion densities and velocities
     &             NHEAT(1:CTIPDIM), & !.. OUT: array, Neutral heating rate (eV/cm^3/s) 
     &             EFLAG)  !.. OUT: 2D array, Error Flags

!dbg20110802: 3D multiple-lp run
!if ( sw_debug )  sw_debug=.false.




! output
      DO ipts=1,CTIPDIM
!dbg20120501
         DO jth=1,ISPEC
            plasma_3d(jth,ipts+IN-1,mp) = XIONNX(jth,ipts)
         END DO !jth

!te
         plasma_3d(ISPEC+1,ipts+IN-1,mp) = TE_TIX(3,ipts)
!ti
         DO jth=1,2
            plasma_3d(jth+ISPEC+1,ipts+IN-1,mp) = TE_TIX(jth,ipts)
         END DO !jth
!vi
         DO jth=1,ISPEV
           plasma_3d(jth+ISPEC+3,ipts+IN-1,mp) = XIONVX(jth,ipts)
         END DO !jth
!dbg20120501         plasma_3d(mp,lp)%N_m3( 1:ISPEC,ipts) = XIONNX(1:ISPEC,ipts)
!dbg20120501         plasma_3d(mp,lp)%Te_k(         ipts) = TE_TIX(3      ,ipts)
!dbg20120501         plasma_3d(mp,lp)%Ti_k(       2,ipts) = TE_TIX(2      ,ipts)
!dbg20120501         plasma_3d(mp,lp)%Ti_k(       1,ipts) = TE_TIX(1      ,ipts)
!dbg20120501         plasma_3d4n(ipts+IN-1,mp)%V_ms1( 1:ISPEV     ) = XIONVX(1:ISPEV,ipts)
!dbg20110927
!dbg20110927         plasma_3d(mp,lp)%V_ms1(1:ISPEV,ipts) = XIONVX(1:ISPEV,ipts)
!dbg20110927         plasma_3d(mp,lp)%heating_rate_e_cgs(  ipts)=EHTX(3   ,ipts) 
!dbg20110927         plasma_3d(mp,lp)%heating_rate_i_cgs(2,ipts)=EHTX(2   ,ipts)
!dbg20110927         plasma_3d(mp,lp)%heating_rate_i_cgs(1,ipts)=EHTX(1   ,ipts)
!dbg20110923         IF ( INNO<0 ) &  !when flip calculates NO
!dbg20110923        &  plasma_3d(mp,lp)%NO_m3(      ipts) =   NNOX(        ipts)
      END DO       !DO ipts=1,CTIPDIM
!dbg20110927      NHEAT_cgs(IN:IS,mp) =  NHEAT(        1:CTIPDIM) 

      PRUNIT_dum = PRUNIT
      CALL WRITE_EFLAG(PRUNIT_dum, &  !.. Unit number to print results
     &                  EFLAG)    !.. Error flag array


!nm20110822: no more allocatable arrays
!      IF ( ALLOCATED( ZX  ) )  DEALLOCATE(  ZX &
!     & ,STAT=stat_alloc         )
!      IF ( ALLOCATED( SLX  ) )  DEALLOCATE(  SLX ) !
!      IF ( ALLOCATED( GLX ) )  DEALLOCATE( GLX )
!      IF ( ALLOCATED( BMX  ) )  DEALLOCATE(  BMX ) !
!      IF ( ALLOCATED( GRX  ) )  DEALLOCATE(  GRX ) !
!      IF ( ALLOCATED( OX  ) )  DEALLOCATE(  OX ) !
!      IF ( ALLOCATED( HX  ) )  DEALLOCATE(  HX ) !
!      IF ( ALLOCATED( N2X  ) )  DEALLOCATE(  N2X ) !
!      IF ( ALLOCATED( O2X  ) )  DEALLOCATE(  O2X ) !
!      IF ( ALLOCATED( HEX  ) )  DEALLOCATE(  HEX ) !
!      IF ( ALLOCATED( N4SX  ) )  DEALLOCATE(  N4SX ) !
!      IF ( ALLOCATED( TNX  ) )  DEALLOCATE(  TNX ) !
!      IF ( ALLOCATED( TINFX  ) )  DEALLOCATE(  TINFX ) !
!      IF ( ALLOCATED( UNX ) )  DEALLOCATE( UNX )
!      IF ( ALLOCATED( SZA_dum  ) )  DEALLOCATE( SZA_dum  )
!      IF ( ALLOCATED( NNOX ) )  DEALLOCATE( NNOX )
!      IF ( ALLOCATED( EHTX ) )  DEALLOCATE( EHTX )
!      IF ( ALLOCATED( TE_TIX ) )  DEALLOCATE( TE_TIX )
!      IF ( ALLOCATED( XIONNX ) )  DEALLOCATE( XIONNX )
!      IF ( ALLOCATED( XIONVX ) )  DEALLOCATE( XIONVX )
!      IF ( ALLOCATED( NHEAT  ) )  DEALLOCATE( NHEAT  )
!if ( stat_alloc/=0 ) then
!  print *, ALLOCATED( ZX )
!  print *,"!STOP! ZX DEALLOCATION FAILD! in flux_tube_solver:",stat_alloc,mp,lp,in,is,jminx,jmaxx,ctipdim
!  STOP
!endif

      NULLIFY(IN,IS)
      END SUBROUTINE flux_tube_solver
