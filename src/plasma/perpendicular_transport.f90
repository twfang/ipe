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
      SUBROUTINE perpendicular_transport ( utime, mp,lp )
      USE module_precision
      USE module_input_parameters,ONLY: sw_debug
      IMPLICIT NONE
!--- INPUT ---
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      INTEGER (KIND=int_prec),INTENT(IN) :: mp
      INTEGER (KIND=int_prec),INTENT(IN) :: lp
!---

      REAL(KIND=real_prec) :: phi_t0(2) !magnetic longitude,phi[rad] at T0(previous time step)
      REAL(KIND=real_prec) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
      INTEGER (KIND=int_prec),DIMENSION(2,2) :: mp_t0,lp_t0 !1st dim:ihem;2nd dim:i0/i1
      REAL(KIND=real_prec) :: r0_apex ![meter]

!---

!print "(' mp=',I4,' lp=',I4,3F10.4)", mp,lp,mlon_rad(mp)&
!&,plasma_grid_3d(IN,mp)%GL, plasma_grid_3d(IS,mp)%GL       

! calculate where the flux tube is coming from (semi-lagulangian issue)
!      CALL stepback_mag ( mp,lp &
!!d     &, mlon_rad(mp) &
!!d     &, plasma_grid_3d(IN,mp)%GL, plasma_grid_3d(IS,mp)%GL &
!     &, phi_t0      , theta_t0 ) 
!if(sw_debug) print *,'stepback_mag finished!'

      CALL stepback_mag_R (utime, mp,lp, phi_t0 , theta_t0, r0_apex )
if(sw_debug) print *,'stepback_magR finished!'

!      CALL find_neighbor_grid ( mp,lp  &
!     &, phi_t0  , theta_t0 &
!     &,  mp_t0  ,    lp_t0)
!if(sw_debug) print *,'find_neighbor_grid finished!'


      CALL find_neighbor_grid_R ( mp,lp, phi_t0, theta_t0, r0_apex &
     &, mp_t0,lp_t0 )
if(sw_debug) print *,'find_neighbor_grid R finished!'

! prepare all the parameters along the flux tube by interpolation, in addition to the adiabatic term, compressional term
      CALL interpolate_flux_tube ( mp,lp, phi_t0,theta_t0, r0_apex &
     &, mp_t0,lp_t0 )
if(sw_debug) print *,'interpolate_flux_tube finished!' 

      END SUBROUTINE perpendicular_transport

      SUBROUTINE find_neighbor_grid ( mp,lp &
&, phi_t0 , theta_t0 &
&,  mp_t0 ,    lp_t0 )
      USE module_precision
      USE module_physical_constants,ONLY: rtd
      USE module_FIELD_LINE_GRID_MKS,ONLY:plasma_grid_GL,JMIN_IN,JMAX_IS,mlon_rad,dlonm90km
      USE module_IPE_dimension,ONLY: NMP,NLP
      USE module_input_parameters,ONLY:sw_perp_transport,sw_debug,HaloSize
     IMPLICIT NONE
!--- INPUT ---
      INTEGER (KIND=int_prec),INTENT(IN) :: mp
      INTEGER (KIND=int_prec),INTENT(IN) :: lp
      REAL(KIND=real_prec),INTENT(IN) :: phi_t0(2) !magnetic longitude,phi[rad] at T0(previous time step)
      REAL(KIND=real_prec),INTENT(IN) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
!---local
      INTEGER (KIND=int_prec) :: ihem,ihem_max
      INTEGER (KIND=int_prec),DIMENSION(2,2), INTENT(OUT) :: mp_t0,lp_t0
      INTEGER (KIND=int_prec) :: lp_min,l,mpx,mpp,mpm
!---

!3:both THETA&PHI:transport included, NH/SH flux tubes are moving separately with different ExB drift
!dbg20120509: IF ( sw_perp_transport(mp)==3 ) THEN 
IF ( sw_perp_transport==3 ) THEN 
  ihem_max=2
ELSE
  ihem_max=1
END IF

!NH only for debug purpose
!assume for the moment that flux tubes in NH/SH are moving together....
!if NH flux tube is moving differently from SH flux, run the loop upto ihem=2, then flux tube interpolation should be done separately between NHv.s.SH
which_hemisphere: DO ihem=1,1  !ihem_max
!!!dbg20120125: why mp=19/20 are not working??? 
!!!dbg20120125:  mlon_deg = phi_t0(ihem)*rtd
!!!dbg20120125:  mp_t0(ihem,1) = INT( (mlon_deg/dlonm90km) , int_prec )+1
!!!dbg20120125:  mp_t0(ihem,2) = mp_t0(ihem,1)+1
  mpx_loop: DO mpx=0,NMP
    if(mpx+1 > HaloSize) then
      print*,'mpx+1 > HaloSize in find_neighbor_grid_R',mpx,HaloSize
      print*,'Increase the halo size or take smaller time steps.'
      print*,'Stopping in find_neighbor_grid_R'
    endif
    mpp=mp+mpx
    if(mpp > NMP) mpp= mpp-NMP
    mpm=mp-mpx
    if(mpm < 1) mpm= NMP+mpm
    IF ( mlon_rad(mpp)<=phi_t0(ihem) .AND. phi_t0(ihem)<mlon_rad(mpp+1) ) THEN
      mp_t0(ihem,1) =mpp
      mp_t0(ihem,2) =mpp+1
      EXIT mpx_loop
    END IF
    IF ( mpm>1.and.mlon_rad(mpm)<=phi_t0(ihem) .AND. phi_t0(ihem)<mlon_rad(mpm-1) ) THEN
      mp_t0(ihem,1) =mpm-1
      mp_t0(ihem,2) =mpm
      EXIT mpx_loop
    END IF
  END DO mpx_loop !: DO mpx=0,NMP
!dbg20120125:
if(sw_debug) print *,'dbg20120125! sub-find_neighbor_grid:', mp_t0(ihem,1:2),phi_t0(ihem)*rtd, mlon_rad(mp_t0(ihem,1:2))*rtd
!STOP

!find  lp0_t0:NH
IF (ihem==1) THEN

!check pole regions!
  IF ( theta_t0(ihem) < plasma_grid_GL( JMIN_IN(1),1 ) ) THEN 
   lp_t0(ihem,1)=-999
   lp_t0(ihem,2)=1
   print *,'sub-Fi: mp',mp,' lp',lp,'needs special pole interpolation'
   RETURN
  ELSE IF ( plasma_grid_GL( JMIN_IN(lp),lp ) <= theta_t0(ihem) ) THEN 
    lp_min =lp
  ELSE  !plasma_grid_3d(IN,lp)%GL > theta_t0(ihem) ) THEN 

if(sw_debug) print *,'sub-Fi: check GL NH[deg]',(90.-plasma_grid_GL( JMIN_IN(lp),lp )*rtd)

    lp_min =lp-5 !not sure if 10 is enough???
    if (lp_min<=0 ) lp_min=1
    IF ( plasma_grid_GL( JMIN_IN(lp_min),lp_min ) > theta_t0(ihem) ) THEN 
      lp_min=lp-5
      if (lp_min<=0 ) lp_min=1
      IF ( plasma_grid_GL( JMIN_IN(lp_min),lp_min ) > theta_t0(ihem) ) THEN 
        print *,'sub-Fi:NH !STOP! not sure if this is working???'
        STOP
      END IF
    END IF
  END IF! ( plasma_grid_3d(IN,lp)%GL <= theta_t0(ihem) ) THEN 

!SMS$PARALLEL(dh, l) BEGIN
lp_loop: DO l=lp_min,NLP-1  !nearest point-->EQ
IF ( plasma_grid_GL( JMIN_IN(l),l )<=theta_t0(ihem) .AND. theta_t0(ihem)<plasma_grid_GL( JMIN_IN(l+1),l+1 )  ) THEN
  lp_t0(ihem,1)=l
  lp_t0(ihem,2)=l+1
  EXIT lp_loop
ELSE
  if (l==NLP-1) then
  print *,'sub-Fi:NH: !STOP! could not find the lp',plasma_grid_GL( JMIN_IN(l),l ),theta_t0(ihem)
  STOP
  end if
END IF

END DO lp_loop!: DO i=lp_min,NLP  !nearest point-->EQ
!SMS$PARALLEL END
END IF !(ihem==1) THEN

if(sw_debug) print *,'sub-Fi: mlon', mlon_rad(mp_t0(ihem,1))*rtd, phi_t0(ihem)*rtd, mlon_rad(mp_t0(ihem,2))*rtd, mp_t0(ihem,1:2)
if(sw_debug) print *,'sub-Fi: mlat',plasma_grid_GL( JMIN_IN(lp_t0(ihem,1)),lp_t0(ihem,1) )*rtd, theta_t0(ihem)*rtd, plasma_grid_GL( JMIN_IN(lp_t0(ihem,2)),lp_t0(ihem,2) )*rtd

END DO which_hemisphere!:  DO ihem=1,ihem_max


      END SUBROUTINE find_neighbor_grid

!20111005
! using r0_apex as in GIP
      SUBROUTINE find_neighbor_grid_R ( mp,lp &
&, phi_t0 , theta_t0 &
&, r0_apex &
&,  mp_t0 ,    lp_t0 )
      USE module_precision
      USE module_physical_constants,ONLY: rtd,earth_radius
      USE module_FIELD_LINE_GRID_MKS,ONLY:plasma_grid_GL,JMIN_IN,JMAX_IS,mlon_rad,dlonm90km, plasma_grid_Z
      USE module_IPE_dimension,ONLY: NMP,NLP
      USE module_input_parameters,ONLY:sw_perp_transport,sw_debug,HaloSize
     IMPLICIT NONE
!--- INPUT ---
      INTEGER (KIND=int_prec),INTENT(IN) :: mp
      INTEGER (KIND=int_prec),INTENT(IN) :: lp
      REAL(KIND=real_prec),INTENT(IN) :: phi_t0(2) !magnetic longitude,phi[rad] at T0(previous time step)
      REAL(KIND=real_prec),INTENT(IN) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
      REAL(KIND=real_prec), INTENT(IN) :: r0_apex ![meter]
!---local
      INTEGER (KIND=int_prec) :: ihem,ihem_max
      REAL(KIND=real_prec) :: Z_t0
      INTEGER (KIND=int_prec),DIMENSION(2,2), INTENT(OUT) :: mp_t0,lp_t0
      INTEGER (KIND=int_prec) :: lp_min,l
!SMS$DISTRIBUTE(dh,NLP) BEGIN
      INTEGER (KIND=int_prec) :: midpoint(NLP)
!SMS$DISTRIBUTE END
      INTEGER (KIND=int_prec) :: lp1,lp2,midpoint1,midpoint2,mpx,mpp,mpm
!---

!3:both THETA&PHI:transport included, NH/SH flux tubes are moving separately with different ExB drift
!dbg20120509 IF ( sw_perp_transport(mp)==3 ) THEN 
IF ( sw_perp_transport==3 ) THEN 
  ihem_max=2
ELSE
  ihem_max=1
END IF

!NH only for debug purpose
!assume for the moment that flux tubes in NH/SH are moving together....
!if NH flux tube is moving differently from SH flux, run the loop upto ihem=2, then flux tube interpolation should be done separately between NHv.s.SH
which_hemisphere: DO ihem=1,1  !ihem_max
!!!dbg20120125: why mp=19/20 are not working??? 
!!!dbg20120125:  mlon_deg = phi_t0(ihem)*rtd
!!!dbg20120125:  mp_t0(ihem,1) = INT( (mlon_deg/dlonm90km) , int_prec )+1
!!!dbg20120125:  mp_t0(ihem,2) = mp_t0(ihem,1)+1
  mpx_loop: DO mpx=0,NMP
    if(mpx+1 > HaloSize) then
      print*,'mpx+1 > HaloSize in find_neighbor_grid_R',mpx,HaloSize
      print*,'Increase the halo size or take smaller time steps.'
      print*,'Stopping in find_neighbor_grid_R'
    endif
    mpp=mp+mpx
    if(mpp > NMP) mpp= mpp-NMP
    mpm=mp-mpx
    if(mpm < 1) mpm= NMP+mpm
    IF ( mlon_rad(mpp)<=phi_t0(ihem) .AND. phi_t0(ihem)<mlon_rad(mpp+1) ) THEN
      mp_t0(ihem,1) =mpp
      mp_t0(ihem,2) =mpp+1
      EXIT mpx_loop
    END IF
    IF ( mpm>1.and.mlon_rad(mpm)<=phi_t0(ihem) .AND. phi_t0(ihem)<mlon_rad(mpm-1) ) THEN
      mp_t0(ihem,1) =mpm-1
      mp_t0(ihem,2) =mpm
      EXIT mpx_loop
    END IF
  END DO mpx_loop !: DO mpx=1,NMP
!dbg20120125:
if(sw_debug) print *,'dbg20120125! sub-find_neighbor_grid_R:', mp_t0(ihem,1:2),phi_t0(ihem)*rtd, mlon_rad(mp_t0(ihem,1:2))*rtd
!STOP

!find  lp0_t0:NH
IF (ihem==1) THEN

!check pole regions!
!not totally sure whether I should use theta_t0 or r0_apex???
IF ( theta_t0(ihem) < plasma_grid_GL( JMIN_IN(1),1 ) ) THEN 
   lp_t0(ihem,1)=-999
   lp_t0(ihem,2)=1
   print *,'sub-Fi: mp',mp,' lp',lp,'needs special pole interpolation'
   RETURN
ELSE IF ( theta_t0(ihem) > plasma_grid_GL( JMIN_IN(NLP),NLP ) ) THEN
   print *,'sub-Fi: !STOP! invalid theta_t0',mp,lp,theta_t0(ihem),plasma_grid_GL( JMIN_IN(NLP),NLP )
   STOP
ELSE   !IF ( plasma_grid_GL( JMIN_IN(lp),lp ) <= theta_t0(ihem) ) THEN 

!!!UNDERCONSTRUCTION!!!
if(sw_debug) print *,'sub-FiR: check GL NH[deg]',(90.-plasma_grid_GL( JMIN_IN(lp),lp )*rtd)

!    lp_min =lp-5 !not sure if 10 is enough???
!    if (lp_min<=0 ) lp_min=1
!    IF ( plasma_grid_GL( JMIN_IN(lp_min),lp_min ) > theta_t0(ihem) ) THEN 
!      lp_min=lp-5
!      if (lp_min<=0 ) lp_min=1
!      IF ( plasma_grid_GL( JMIN_IN(lp_min),lp_min ) > theta_t0(ihem) ) THEN 
!        print *,'sub-Fi:NH !STOP! not sure if this is working???'
!        STOP
!      END IF
!    END IF

!SMS$PARALLEL(dh, l) BEGIN
lp_loop0: DO l=1,NLP  !longest -->shortest flux tube
 midpoint(l) = JMIN_IN(l) + ( JMAX_IS(l) - JMIN_IN(l) )/2
END DO lp_loop0

z_t0 = r0_apex - earth_radius


!d l=130
!d print *,JMIN_IN(l),JMAX_IS(l), midpoint(l),z_t0


lp_loop: DO l=1,NLP-1  !longest -->shortest flux tube
 
  IF ( plasma_grid_Z( midpoint(l+1),l+1 )<=Z_t0 .AND. Z_t0<plasma_grid_Z( midpoint(l),l )  ) THEN
    lp_t0(ihem,1)=l   !1outer flux tube
    lp_t0(ihem,2)=l+1 !2inner flux tube
    EXIT lp_loop
  ELSE
    if (l==NLP-1) then
      print *,'sub-FiR:NH: !STOP! could not find the lp',plasma_grid_Z( midpoint(l),l ),z_t0
      STOP
    end if
  END IF

END DO lp_loop!: DO i=lp_min,NLP  !nearest point-->EQ
!SMS$PARALLEL END

!OUT  lp_t0(ihem,1)=l  
!IN   lp_t0(ihem,2)=l+1
lp1 = lp_t0(ihem,1)
midpoint1 = midpoint(lp1)
lp2 = lp_t0(ihem,2) !=l+1
midpoint2 = midpoint(lp2)
if(sw_debug) print *,lp1,midpoint1,lp2,midpoint2
if(sw_debug) print *,'sub=FiR: check R',(earth_radius+plasma_grid_Z(midpoint1,lp1)),(earth_radius+Z_t0),(earth_radius+plasma_grid_Z(midpoint2,lp2))

if(sw_debug) print *,'sub-FiR: mlon',mlon_rad(mp_t0(ihem,1))*rtd, phi_t0(ihem)*rtd, mlon_rad(mp_t0(ihem,2))*rtd, mp_t0(ihem,1:2)
if(sw_debug) print *,'sub-FiR: mlat',(90.-plasma_grid_GL( JMIN_IN(lp_t0(ihem,1)),lp_t0(ihem,1) )*rtd), (90.-theta_t0(ihem)*rtd),(90.- plasma_grid_GL( JMIN_IN(lp_t0(ihem,2)),lp_t0(ihem,2) )*rtd)


END IF! ( plasma_grid_3d(IN,lp)%GL <= theta_t0(ihem) ) THEN 
END IF !(ihem==1) THEN
END DO which_hemisphere!:  DO ihem=1,ihem_max

if(sw_debug) print *,  'sub-find_neighbor_grid R finished!'
      END SUBROUTINE find_neighbor_grid_R
