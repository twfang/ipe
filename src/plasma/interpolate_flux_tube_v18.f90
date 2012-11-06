!dbg20120330: next version! should use the apex routine to precisely estimate x(0),r(0), b0(0) for the imaginary flux tube...(one change at a time...)
!           :also density interpolation should be done with log???
!v17: previous versions were all WRONG!!! regarding how to implement ksi factor...
!v11: 20111101: ksi fac 1.0
!v15: 20111212: included Te/i transport
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
      SUBROUTINE interpolate_flux_tube (mp,lp &
&, phi_t0,theta_t0 &
&, r0_apex &
&,mp_t0,lp_t0 )
      USE module_precision
!     plasma_grid_3d,plasma_grid_Z,plasma_grid_GL,plasma_3d_old are all IN arrays
      USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,JMAX_IS,plasma_grid_3d,plasma_grid_Z,plasma_grid_GL,ht90,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,plasma_3d_old
      USE module_input_parameters,ONLY:sw_perp_transport,sw_debug,sw_ksi
      USE module_plasma,ONLY:plasma_1d !dbg20120501 n0_1d   !d, n0_2dbg
      USE module_IPE_dimension,ONLY: ISPEC,ISPET,IPDIM
      USE module_physical_constants,ONLY: earth_radius,pi,zero
      IMPLICIT NONE
!--- INPUT ---
      INTEGER (KIND=int_prec),INTENT(IN) :: mp
      INTEGER (KIND=int_prec),INTENT(IN) :: lp
      REAL(KIND=real_prec),INTENT(IN) :: phi_t0(2) !magnetic longitude,phi[rad] at T0(previous time step)
      REAL(KIND=real_prec),INTENT(IN) :: theta_t0(2) !magnetic latitude,theta[rad] at T0
      INTEGER (KIND=int_prec),DIMENSION(2,2), INTENT(IN) :: mp_t0,lp_t0 !1st dim:ihem, 2nd dim:ilp/imp
      REAL(KIND=real_prec),INTENT(IN) :: r0_apex ![meter]
!---local
      INTEGER (KIND=int_prec) :: imp_max,ihem_max, ihem
      INTEGER (KIND=int_prec) :: ilp,lp0,imp,mp0
      INTEGER (KIND=int_prec) :: ip,i,jth
      INTEGER (KIND=int_prec) :: ispecial,isouth,inorth,i1d,ip1d,midpoint ,kk
      REAL(KIND=real_prec) :: factor2

      REAL(KIND=real_prec) :: factor
      REAL(KIND=real_prec), DIMENSION(2) :: ksi_fac  !dim:ilp
      REAL(KIND=real_prec),DIMENSION(0:2) :: r,lambda_m,rapex,B0,x
      INTEGER (KIND=int_prec),PARAMETER :: TSP=3      !N(1:3) perp.transport
      INTEGER (KIND=int_prec),PARAMETER :: iT=ISPEC+3   !add T(1:3)
      REAL(KIND=real_prec) :: n0(iT,2) !1d:species(N&T); 2dim:ilp

      INTEGER (KIND=int_prec),PARAMETER :: iB=iT+1 !add B
      INTEGER (KIND=int_prec),PARAMETER :: iR=iB+1 !add R
      REAL(KIND=real_prec) :: Qint(iR,IPDIM,2,2)  !1d:species; ; 4d:imp; 5d:ilp
!---

!array initialization: may not be necessary because they are local parameters...
      Qint(:,:,:,:)=zero
      n0(:,:)=zero

!dbg
if(sw_debug) print *,'lp_t0',lp_t0
if(sw_debug) print *,'mp_t0',mp_t0

!dbg20120509: IF ( sw_perp_transport(mp)==1 ) THEN !THETA only transport included
IF ( sw_perp_transport==1 ) THEN !THETA only transport included
  imp_max=1  
!dbg20120509 ELSE IF ( sw_perp_transport(mp)>=2 ) THEN
ELSE IF ( sw_perp_transport>=2 ) THEN
  imp_max=2
END IF

!3:both THETA&PHI:transport included, NH/SH flux tubes are moving separately with different ExB drift
!dbg20120509: IF ( sw_perp_transport(mp)==3 ) THEN 
IF ( sw_perp_transport==3 ) THEN 
  ihem_max=2
ELSE
  ihem_max=1
END IF



!NH only for debug purpose
!let us assume for the moment that flux tubes in NH/SH are moving together with the same ExB drift....
!if NH flux tube is moving differently from SH flux, run the loop upto ihem=2, then flux tube interpolation should be made separately between NHv.s.SH
which_hemisphere: DO ihem=1,1 !ihem_max
if(sw_debug) print *,'ihem',ihem

!check if the flux tube needs special pole interpolation???
IF ( lp_t0(ihem,1)>=1 ) THEN

 

! grid point distributions are identical in all mp: between mp_t0(ihem,1)v.s.mp_t0(ihem,2)
! however Q values are not identical in all mp!!! 
mp_t0_loop: DO imp=1,imp_max
  mp0 = mp_t0(ihem,imp)
if(sw_debug) print *,ihem,imp,'mp_t0',mp_t0(ihem,imp),' mp0',mp0

!(1) Q interpolation: from Q_T0(mp_t0,lp_t0) --> Q_T1(mp,lp) for all 4 flux tubes
  lp_t0_loop: DO ilp=1,2 !outer/inner flux tubes
if(sw_debug) print *,'ilp',ilp
if(sw_debug) print *,'lp_t0',lp_t0(ihem,ilp)

! (1.1) i0:  lp_t0(ihem,1)=l
    lp0 = lp_t0(ihem,ilp)
if(sw_debug) print *,'lp0',lp0
if(sw_debug) print *,'mp',mp,'lp',lp

    flux_tube_loopT1_Q: DO ip=JMIN_IN(lp),JMAX_IS(lp)
      ip1d=ip-JMIN_IN(lp)+1
if(sw_debug)  print *,'ip=',ip,' ip1d=',ip1d
!check the foot point values
!ispecial=2: interpolation at/below IN

! Will this be in the halo? JFM
      IF ( plasma_grid_3d( JMIN_IN(lp0) , lp0,mp0,IQ) < plasma_grid_3d(ip ,lp,mp,IQ) ) THEN
       ispecial=2
       isouth=JMIN_IN(lp0)+1  !not used!!!
       inorth=JMIN_IN(lp0)
       i1d  =+1


if ( sw_debug.and.lp0==149 .and. ip1d>=63 ) then 
print *,'!dbg20120508! ispecial=',ispecial,ip,i,JMAX_IS(lp0),lp0,mp0, plasma_grid_3d(i,lp0,mp0,IQ) ,plasma_grid_3d(ip,lp,mp,IQ),lp,mp
endif


!ispecial=1: interpolation at/below IS
      ELSE IF ( plasma_grid_3d( JMAX_IS(lp0) , lp0,mp0,IQ) >= plasma_grid_3d(ip,lp,mp,IQ) ) THEN
       ispecial=1
       isouth=JMAX_IS(lp0)   !not used!!!
       inorth=JMAX_IS(lp0)-1
       i1d  =JMAX_IS(lp0)-JMIN_IN(lp0)+1



if ( sw_debug.and.lp0==149 .and. ip1d>=63 ) then 
print *,'!dbg20120508! ispecial=',ispecial,ip,i,JMAX_IS(lp0),lp0,mp0, plasma_grid_3d(i,lp0,mp0,IQ) ,plasma_grid_3d(ip,lp,mp,IQ),lp,mp
endif


!dbg20120504!!???why this output never appear on ipeXXX.log??? 
!???why there is never ispecial=1???
if(sw_debug.and.lp==149)then
print *,'!dbg20120504 Q=',plasma_grid_3d( JMAX_IS(lp0) , lp0,mp0,IQ) , plasma_grid_3d(ip,lp,mp,IQ) ,JMAX_IS(lp0) , lp0,mp0, ip,ip1d,lp,mp,isouth,inorth,i1d
endif

      ELSE

!search for north & south grid point of i1
        flux_tube_loopT0: DO i=JMIN_IN(lp0),JMAX_IS(lp0)

         IF ( plasma_grid_3d(i,lp0,mp0,IQ) < plasma_grid_3d(ip,lp,mp,IQ) ) THEN
! ispecial=0: normal interpolation
           ispecial = 0
           inorth = i-1
           isouth = i
           i1d   = i-JMIN_IN(lp0)+1
         
! factor2 IS NOT equal for all mp
           factor2=(plasma_grid_3d(ip,lp,mp,IQ)-plasma_grid_3d(isouth,lp0,mp0,IQ))/(plasma_grid_3d(inorth,lp0,mp0,IQ)-plasma_grid_3d(isouth,lp0,mp0,IQ))
if ( factor2<0.0.or.factor2>1.0) then
WRITE(6,*) 'sub-Intrp:!STOP! invalid factor2',factor2 ,ip,lp,mp,isouth,lp0,mp0
STOP
endif       



if ( sw_debug.and.lp0==149 .and. ip1d>=63 ) then 
print *,'!dbg20120508!',ip,factor2,i,JMAX_IS(lp0),lp0,mp0, plasma_grid_3d(i,lp0,mp0,IQ) ,plasma_grid_3d(ip,lp,mp,IQ),lp,mp
endif

           EXIT flux_tube_loopT0
         END IF 




        END DO flux_tube_loopT0 !: DO i=IN,IS

      END IF !( Q_t0(IN) < Q_t1(ip) ) THEN


!ispecial=0: normal interpolation
      if(ispecial == 0) then
! calculate all the ionospheric parameters      
!        ni1_in(ip)=(factor2*(ni(inorth,mp0,1) - ni(isouth,mp0,1))) + ni(isouth,mp0,1)


!N 1:TSP: density
!not sure if LOG is necessary for densities???
jth_loop0:         DO jth=1,iT !=TSP+3

IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop0
!dbg20120501            IF(jth<=TSP) THEN
               Qint(jth, ip1d,imp,ilp) = (factor2*(plasma_3d_old(inorth,lp0,mp0,jth) - plasma_3d_old(isouth,lp0,mp0,jth))) + plasma_3d_old(isouth,lp0,mp0,jth)

!T TSP+1:TSP+3=iT
!dbg20120501            ELSE IF(jth==TSP+1) THEN
!dbg20120501               Qint(jth, ip1d,imp,ilp) = (factor2*(plasma_3d_old(mp0,lp0)%Te_k(i1d-1) - plasma_3d_old(mp0,lp0)%Te_k(i1d))) + plasma_3d_old(mp0,lp0)%Te_k(i1d)

!dbg20120501            ELSE !Ti
!dbg20120501               Qint(jth, ip1d,imp,ilp) = (factor2*(plasma_3d_old(mp0,lp0)%Ti_k( (jth-TSP-1),i1d-1) - plasma_3d_old(mp0,lp0)%Ti_k( (jth-TSP-1),i1d))) + plasma_3d_old(mp0,lp0)%Ti_k( (jth-TSP-1),i1d)
!dbg20120501            END IF

            if (&
!&jth==1&
&jth==5&
&.and.Qint(jth, ip1d,imp,ilp)<=0.) then
!dbg20120501            if (jth==1.and.Qint(jth, ip1d,imp,ilp)<=0.) then



               WRITE(6,*)'sub-Intrp:!STOP! INVALID density',Qint(jth, ip1d,imp,ilp),factor2 &
                    &,plasma_3d_old(inorth,lp0,mp0,jth)   & !dbg20120501
                    &,plasma_3d_old(isouth,lp0,mp0,jth)   & !dbg20120501
                    &,jth, ip1d,imp,ilp,mp0,lp0,i1d,inorth,isouth
               STOP
            endif
         END DO jth_loop0 !jth=1,iT !=TSP+3

!B iT+1=iB: B magnetic field intensity
         Qint(iB, ip1d,imp,ilp) = (factor2*(plasma_grid_3d(inorth,lp0,mp0,IBM) - plasma_grid_3d(isouth,lp0,mp0,IBM))) + plasma_grid_3d(isouth,lp0,mp0,IBM)

!NOTE: R should be the same for all mp!!!
!R iB+1=iR: R = RE + Z
         Qint(iR, ip1d,imp,ilp) =( (factor2*(plasma_grid_Z(inorth,lp0) - plasma_grid_Z(isouth,lp0))) + plasma_grid_Z(isouth,lp0) ) +earth_radius

if ( sw_debug.and.lp0==149.and.ip1d>=63) then
print *,'!dbg20120504 R=' &
, Qint(iR, ip1d,imp,ilp) &
, (factor2*(plasma_grid_Z(inorth,lp0) - plasma_grid_Z(isouth,lp0))) &
&, plasma_grid_Z(inorth,lp0) &
&, plasma_grid_Z(isouth,lp0) &
&, factor2 &
&,lp0,ip1d,inorth,isouth

endif






!ispecial=1: interpolation at/below IS_t0
      ELSE if(ispecial == 1) then

        jth_loop1: DO jth=1,iT
IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop1
          Qint(jth   ,ip1d,imp,ilp) = plasma_3d_old(isouth,lp0,mp0,jth)
        END DO jth_loop1 !jth
         !N:        ni1_in(ip)=ni(IS_t0,mp0,1)
!dbg20120501        Qint(1:TSP   ,ip1d,imp,ilp) = plasma_3d_old(mp0,lp0)%N_m3( 1:TSP,i1d)
        !Te:
!dbg20120501         Qint(TSP+1   ,ip1d,imp,ilp) = plasma_3d_old(mp0,lp0)%Te_k(         i1d)
        !Ti:
!dbg20120501         Qint(TSP+2:iT,ip1d,imp,ilp) = plasma_3d_old(mp0,lp0)%Ti_k(1:ISPET,i1d)
        !B:
        Qint(iB      ,ip1d,imp,ilp) = plasma_grid_3d(isouth,lp0,mp0,IBM)
        !R:
        Qint(iR      ,ip1d,imp,ilp) = plasma_grid_Z(isouth,lp0) +earth_radius

!dbg20120504!!???why this output never appear on ipeXXX.log??? 
if(sw_debug.and.lp==149) then
print *,'!dbg20120504',        Qint(iR      ,ip1d,imp,ilp),plasma_grid_Z(isouth,lp0) ,isouth,JMIN_IN(lp0),JMAX_IS(lp0),lp0,mp0,iR,ip1d,imp,ilp
endif

!ispecial=2: interpolation at/below IN_t0
     ELSE if(ispecial == 2) then
       jth_loop2: DO jth=1,iT
IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop2
          Qint(jth   ,ip1d,imp,ilp) = plasma_3d_old(inorth,lp0,mp0,jth)
       END DO jth_loop2!jth
        !N        ni1_in(ip)=ni(IN_t0,mp0,1) 
!dbg20120501        Qint(1:TSP   ,ip1d,imp,ilp) = plasma_3d_old(mp0,lp0)%N_m3(1:TSP,i1d)
         !Te
!dbg20120501        Qint(TSP+1   ,ip1d,imp,ilp) = plasma_3d_old(mp0,lp0)%Te_k(        i1d)
         !Ti
!dbg20120501        Qint(TSP+2:iT,ip1d,imp,ilp) = plasma_3d_old(mp0,lp0)%Ti_k(1:ISPET,i1d)
         !B
        Qint(iB      ,ip1d,imp,ilp) = plasma_grid_3d(inorth,lp0,mp0,IBM)
         !R
        Qint(iR      ,ip1d,imp,ilp) = plasma_grid_Z(inorth,lp0) +earth_radius

      END IF !(ispecial == 1) then
    END DO flux_tube_loopT1_Q !: DO ip=IN,IS
  END DO lp_t0_loop !: DO ilp=1,2
END DO mp_t0_loop!: DO imp=1,2





! (2) intepolate between the 4 flux tubes each point with the same Q value
! what kind of interpolation is the most appropriate for the 4 points???
!note: 2 flux tubes only used if imp_max=1/sw_perp_transport==1 ---THETA only transport included
!I should check quadratic interpolation CTIPe has using the third flux tube: tubes_quadratic_interpolate.f 
!here in this interpolation should be done for both mp=j0,j1 simultaneously...

! compute factor
! APEX latitude[rad]: eq(3.3) of the imaginary FT(phi0,theta0)
!dbg20120503: lambda_m(0) = pi*0.5 - theta_t0(ihem)
rapex(0) = r0_apex
lambda_m(0) = ACOS(SQRT((earth_radius+ht90)/r0_apex  ))

!R of apex altitude for the two IN/OUT FTs
DO ilp=1,2  !outer/inner flux tubes

if(sw_debug) print *,'!dbg20120503 lp=',lp_t0(ihem,ilp),ilp,ihem

  lp0 = lp_t0(ihem,ilp)
  midpoint = JMIN_IN(lp0) + ( JMAX_IS(lp0) - JMIN_IN(lp0) )/2
  rapex(ilp)=plasma_grid_Z(midpoint,lp0) + earth_radius

if(sw_debug) print *,'!dbg20120503! rapex',rapex(ilp),midpoint
  
!note: this factor cannot work when lp<=6!!!
  IF ( lp>6 ) THEN
    lambda_m(ilp)= ACOS(SQRT((earth_radius+ht90)/rapex(ilp)))
  ELSE
  !note: this factor cannot work when lp<=6!!!
  !because rapex does not mean anything for lp<=6
    lambda_m(ilp)  = pi*0.5 - plasma_grid_GL( JMIN_IN(lp0),lp0 )
  END IF
END DO  !DO ilp=1,2
  
!not sure which factor is more correct??? either r- or lambda base???
IF ( lp>6 .or. rapex(1)==rapex(2) ) THEN
! r = RE + ha(APEX height)
!!!  rapex(0)=( earth_radius + ht90 ) / ( COS( lambda_m(0) ) * COS( lambda_m(0)) )

  factor = ( rapex(0)-rapex(2) ) / ( rapex(1)-rapex(2) )
ELSE !IF lp<=6
!???not sure if the rapex(0) mean anything for huge flux tubes??? thus use the factor of the magnetic apex latitude as in GIP...
! the values are only for NH

  factor = ( lambda_m(0) - lambda_m(2)) / (lambda_m(1) - lambda_m(2))
if(sw_debug)  print *,'!!!different factor!!!',factor,mp,lp
END IF

!error trap
IF ( factor>1.0.OR.factor<0.0) THEN
  print *,'!!!INVALID factor!!!',factor,mp,lp,rapex(0:2),lambda_m(0:2)
  IF ( factor>1.0 )factor=1.0
  IF ( factor<0.0 )factor=0.0
ENDIF



!!!CAUTION!!! this mp loop does not work when imp_max=2!!!
mp_t0_loop1: DO imp=1,imp_max
  flux_tube_loopT1_fac: DO i=JMIN_IN(lp),JMAX_IS(lp) !9000
    i1d=i-JMIN_IN(lp)+1
!1 FTin; 2 FTout
    r(1:2)=Qint(iR,i1d,imp,1:2) !R
if(sw_debug) print "(2i8,'r12=',4E12.4)",i,i1d,Qint(iR,i1d,imp,1:2)

!here what is the best way to calculate r(0)of the imaginary FT??? 
!TODO!!! I should plot this r(0) field line to see if it looks reasonable???
! to be more precise, I could use the apex routine to generate this flux tube, but it would become very computationally expensive... 
    r(0) = factor * ( r(1)-r(2) ) + r(2) 
if(sw_debug) print "('r0=',4E12.4)",factor,r(0:2)


!dbg20120503:
if ( sw_debug.and.lp==149 ) then
write(unit=7000,fmt="(i6,i3,3E13.5)") i,i1d,  Qint(iR,i1d,imp,1), r(0),Qint(iR,i1d,imp,2)
write(unit=7001,fmt="(i6,i3,5F11.5)") i,i1d, (Qint(iR,i1d,imp,1)-earth_radius)*1.0E-3, (r(0)-earth_radius)*1.0e-3,(Qint(iR,i1d,imp,2)-earth_radius)*1.0e-3, plasma_grid_z(i,lp)*1.0e-3, plasma_grid_3d(i,lp,mp,IQ) 
endif



!weighting of X between Nin & Nout
! X can be either R or lambda (but only at IN/IS!!!)
IF( r(1)/=r(2) ) then
   x(0:2) = r(0:2)
if(sw_debug) print *,'!dbg20120503: x',x(0:2)
ELSE IF( r(1)==r(2) ) THEN
   if(sw_debug) print "('!R1=R2!',5E13.5,i6)",r(0:2),Qint(iR,i1d,imp,1:2),i
   IF( r(2)==(earth_radius+ht90) ) then
!dbg20111006: somehow IS does not fit here???
     if(sw_debug)  print *,'i=IN/S',i,i1d,JMIN_IN(lp),JMAX_IS(lp),lambda_m(0:2)*180./pi
     x(0:2) = lambda_m(0:2)
if(sw_debug) print *,'!dbg20120503: lambda',x(0:2)
   ELSE
     WRITE(6,*)'sub-Intrp:!STOP! INVALID R12!!!',iR,i1d,imp
     STOP
   END IF
END IF




! 1. interpolate Bfield intensity Bt0 at the imaginary FT(phi0,theta0) using dipole assumption
if(sw_debug) print "('QintB=',3E12.4)",Qint(iB,i1d,imp,1:2)
    B0(1:2)=Qint(iB,i1d,imp,1:2) * ( r(1:2)*r(1:2)*r(1:2) )/(r(0)*r(0)*r(0))
    B0(0)=( B0(1)+B0(2) )*0.50
if(sw_debug) print "('B=',3E12.4)",B0(0),B0(1),B0(2)
!dbg20111101:v14
    B0(0) = ( (x(1)-x(0))*Qint(iB,i1d,imp,2) + (x(0)-x(2))*Qint(iB,i1d,imp,1) ) / ( x(1)-x(2) )
if(sw_debug) print "('v14:B=',3E12.4)",B0(0) !,B0(1),B0(2)


!dbg20120503:
if (sw_debug.and. lp==149 ) then
write(unit=7002,fmt="(i6,i3,4E13.5)") i,i1d, Qint(iB,i1d,imp,1), B0(0),Qint(iB,i1d,imp,2),plasma_grid_3d(i,lp,mp,IBM)
endif





    if ( sw_ksi==0 ) then
      ksi_fac(1:2) =1.000
    else if ( sw_ksi==1 ) then
      ksi_fac(1:2) = B0(0) / Qint(iB,i1d,imp,1:2)
WRITE(6,*)'!STOP! INVALID option! sw_ksi=',sw_ksi
STOP

!dbg20120330: new and CORRECT method to apply the ksi factor!
    else if ( sw_ksi==2 ) then
      ksi_fac(1) = plasma_grid_3d(i,lp,mp,IBM) / B0(0) 
    end if
if(sw_debug) print "('ksi=',2E12.4)", ksi_fac(1:2)*ksi_fac(1:2)

!???1:in; 2:out???
    jth_loop3: DO jth=1,iT !=TSP+3
IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop3
!dbg20120330      IF ( jth<=TSP ) THEN
        !N interpolate from Nin onto FT(phi0,theta0) by applying ksi factor^2 eq(9) p117 PGR thesis
        n0(jth,1:2)=Qint(jth,i1d,imp,1:2) !dbg20120330* ksi_fac(1:2)*ksi_fac(1:2)
!dbg20120330      ELSE !IF ( jth>=TSP+1 ) THEN         
!dbg20120330        !T: ksi factor^4/3: eq (9) page 117 in PGR thesis for T
!dbg20120330        n0(jth,1:2)=Qint(jth,i1d,imp,1:2) * ksi_fac(1:2)**(4./3.)
!dbg20120330      END IF
    END DO jth_loop3!jth=1,TSP+3


!dbg20120501: near future two 
!(1) n0(jth,1:2) can be replaced by Qint(,,,) in below lines and above lines can be commented!
!(2)  IF ( jth<=TSP ) THEN  !for densities
! factor1 = ksi_fac(1)*ksi_fac(1)
! ELSE !ID(jth>TSP) THEN
! factor1 = ksi_fac(1)**(4./3.) 
! END IF
! plasma_1d(jth,i1d) = factor1 * ( (x(1)-x(0))*n0(jth,2) + (x(0)-x(2))*n0(jth,1) ) / ( x(1)-x(2) )


!4. calculate N(phi0,theta0) with weighting of X between Nin & Nout
! X can be either R or lambda (but only at IN/IS)
!    plasma_3d_old(mp,lp)%N_m3(1:TSP,i1d) &
    jth_loop4: DO jth=1,iT !TSP+3
IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop4

       IF ( (x(1)-x(2))/=0.) THEN
!dbg20120501          IF ( jth<=TSP ) THEN
!dbg20120501             n0_1d%N_m3(jth        ,i1d) = ( (x(1)-x(0))*n0(jth,2) + (x(0)-x(2))*n0(jth,1) ) / ( x(1)-x(2) )
             plasma_1d(jth,i1d) = ( (x(1)-x(0))*n0(jth,2) + (x(0)-x(2))*n0(jth,1) ) / ( x(1)-x(2) )

       !error check
!dbg20120501             IF (n0_1d%N_m3(jth,i1d)<=0.) THEN
             IF (plasma_1d(jth,i1d)<=0.) THEN

if(sw_debug) print *,'!dbg20120503' &
&, ( (x(1)-x(0))*n0(jth,2) + (x(0)-x(2))*n0(jth,1) ), ( x(1)-x(2) )&
&, ( (x(1)-x(0))*n0(jth,2)) , ( (x(0)-x(2))*n0(jth,1) ) &
&,  (x(1)-x(0)) , n0(jth,2),  (x(0)-x(2)), n0(jth,1) &
&,   x(0), x(1),x(2) 

                print "('sub-int:!STOP! INVALID density/temp',3E12.4,5i7)" & 
                     &, n0(jth,1), plasma_1d(jth,i1d) , n0(jth,2),jth,i1d,i,mp,lp
                print "('!check X!=',3E12.4)",x(1),x(0),x(2)
                print "('!check B!=',3E12.4)",B0(1),B0(0),B0(2)
                STOP
             END IF

             IF ( jth<=TSP ) THEN
!dbg20120330: new ksi_factor
                plasma_1d(jth,i1d) = plasma_1d(jth,i1d) *(ksi_fac(1)*ksi_fac(1))

!dbg20120501          ELSE IF ( jth==TSP+1 ) THEN 
!dbg20120501             n0_1d%Te_k(            i1d) = ( (x(1)-x(0))*n0(jth,2) + (x(0)-x(2))*n0(jth,1) ) / ( x(1)-x(2) )


       !error check with Te
!dbg20120501             IF (n0_1d%Te_k(i1d)<=0.) THEN
!dbg20120501                print "('sub-int:!STOP! INVALID Te',3E12.4,5i7)" & 
!dbg20120501                     &, n0(jth,1), n0_1d%Te_k(i1d) , n0(jth,2),jth,i1d,i,mp,lp
!dbg20120501                print "('!check X!=',3E12.4)",x(1),x(0),x(2)
!dbg20120501                print "('!check B!=',3E12.4)",B0(1),B0(0),B0(2)
!dbg20120501                STOP
!dbg20120501             END IF

!dbg20120330: new ksi_factor
!dbg20120501 plasma_1d(jth,i1d) = plasma_1d(jth,i1d) *(ksi_fac(1)**(4./3.))

             ELSE !             IF ( jth>TSP ) THEN
!dbg20120501             n0_1d%Ti_k(jth-TSP-1,i1d) = ( (x(1)-x(0))*n0(jth,2) + (x(0)-x(2))*n0(jth,1) ) / ( x(1)-x(2) )

       !error check with Tion
!dbg20120501             IF (n0_1d%Ti_k(jth-TSP-1,i1d)<=0.) THEN
!dbg20120501                print "('sub-int:!STOP! INVALID Ti',3E12.4,5i7)" & 
!dbg20120501                     &, n0(jth,1), n0_1d%Ti_k(jth-TSP-1,i1d) , n0(jth,2),jth,i1d,i,mp,lp
!dbg20120501                print "('!check X!=',3E12.4)",x(1),x(0),x(2)
!dbg20120501                print "('!check B!=',3E12.4)",B0(1),B0(0),B0(2)
!dbg20120501                STOP
!dbg20120501             END IF

!dbg20120330: new ksi_factor
                plasma_1d(jth,i1d) = plasma_1d(jth,i1d) *(ksi_fac(1)**(4./3.))

          END IF !             IF ( jth<=TSP ) THEN
       ELSE
          WRITE(6,*)'sub-Intrp:!STOP! INVALID x(1:2)',x(0:2),i1d,mp,lp
          STOP
       END IF



    END DO jth_loop4!jth=1,iT !=TSP+3



 END DO flux_tube_loopT1_fac !: DO i=in(lp),is(lp) !9000
 
END DO mp_t0_loop1 !: DO imp=1,imp_max


ELSE IF ( lp_t0(ihem,1)==-999 ) THEN
 print *,'sub-intrp:!STOP! sub-Int: mp',mp,' lp',lp,'needs special pole interpolation: right now we do not have such special interpolation...'
 STOP
ELSE
print *,'sub-Intrp:!STOP! INVALID lp_t0:',lp_t0,' mp',mp,' lp',lp
STOP
END IF

      END DO which_hemisphere !: DO ihem=1,ihem_max



      END SUBROUTINE interpolate_flux_tube
