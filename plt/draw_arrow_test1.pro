pro draw_arrow_test1,u,v,thetaR,radi $
,sunlons1,nmp,nlp
;,X00,dX,X_SIZE,Y0,Y_SIZE


ph0=fltarr(nmp,nlp)
ph1=fltarr(nmp,nlp)
th0=fltarr(nmp,nlp)
th1=fltarr(nmp,nlp)
read_fort9001, ph0,ph1,th0,th1




size_resultu=SIZE(u)
print,' size_u',size_resultu
size_resultv=SIZE(v)
print,' size_v',size_resultv
size_radi=SIZE(radi)
print,' size_radi', size_radi
size_thetaR=SIZE(thetaR)
print,' size_thetaR', size_thetaR

size_ph0=SIZE(ph0)
print,' size_ph0', size_ph0








ArrowRef=300. ;2000. ;m/s
;loadct, 0
ArrowCol=0.;255. ;defalt black ;256 ;
factor=0.025 ;0.05
value_thick=0.8

istep=1
jstep=1

for i=0L,size_thetaR[1]-1,   istep  do begin ;thetaR
  for j=1L,  size_radi[1]-1, jstep  do begin ;radi

;dbg20150310
;d if j ne 20 then continue

;ph0-->mlt0
   mlt0    = ph0[i,j]/!DTOR/15.0D0 - sunlons1 * 12.0D0 / !PI   +12.0 ;[hr] !CORRECT!

  if ( mlt0 lt  0. ) then  mlt0 = mlt0 MOD 24.
  if ( mlt0 ge 24. ) then  mlt0 = mlt0 - 24.

  mlt0 = mlt0*!PI/12.0D0   ;MLT_hr --> THETA[rad]

;shift the MLT so that 00MLT comes at the bottom of the plot!
;clockwise 90deg rotation
  mlt0 = mlt0 - !PI*0.50D0       ;(radian)


  ;print, 'check thetaR',i, (thetaR[i]/!DTOR),thetaR[i]

;only NH
 ; if ( which_hem eq 'NH' ) then begin
    comlat0 = $
;90. - 
th0[i,j]/!DTOR ;[deg]
 ; endif else if ( which_hem eq 'SH' ) then  begin
 ;   comlat0 = 90. + th0[i,j]/!DTOR ;[deg] ;degrees measured from South Pole
 ; endif


;if ( radi[j] gt 40. ) then begin
;ArrowCol=255.;0. ;black 
value_thick=2.
;endif else begin
;ArrowCol=255.
;value_thick=0.8
;endelse

X0_arrow = comlat0 * COS(mlt0)
Y0_arrow = comlat0 * SIN(mlt0)

U0=u[i,j]
V0=v[i,j]
mag = SQRT(U0*U0 + V0*V0)
Xmag = V0*cos(thetaR[i])-U0*sin(thetaR[i])
Ymag = V0*sin(thetaR[i])+U0*cos(thetaR[i])
costheta = Xmag/mag
if Ymag ge 0 then $
  sintheta = SQRT(1.0 -costheta*costheta) $
else $
  sintheta = - SQRT(1.0 -costheta*costheta)

comlat_min = -59.  +90. ;lat [deg.] from the pole
comlat_max = -50.  +90.
mlt_rad_min = (270.-15.*3.)*!DTOR
mlt_rad_max = (270.+15.*3.-360.)*!DTOR

;print, 'mlt_rad_min', mlt_rad_min,'mlt_rad_max',mlt_rad_max
;if ( mag gt 800. ) then 
;if ( radi[j] le 10. ) then mag = 0.0
;if ( radi[j] le comlat_min ) or ( radi[j] ge comlat_max ) then CONTINUE
;if ( thetaR[i] le mlt_rad_min ) and ( thetaR[i] ge mlt_rad_max ) then CONTINUE

  arwmag = mag * factor
  dX_arrow = arwmag * costheta
  dY_arrow = arwmag * sintheta

  X1_arrow = radi[j]*COS(thetaR[i]) ;X0_arrow + dX_arrow
  Y1_arrow = radi[j]*SIN(thetaR[i]) ;Y0_arrow + dY_arrow
;print, 'check arrow ref location!!', X0_arrow, Y0_arrow, X1_arrow, Y1_arrow

    ARROW, X0_arrow, Y0_arrow, X1_arrow, Y1_arrow  $
  , /DATA       $ 
;, /NORMALIZED $
  , HSIZE= (!D.X_SIZE / 64.)*0.5 $
  , COLOR=ArrowCol $ ;index] $
;, HTHICK=value] $
;, /SOLID] $
  , THICK=value_thick

  endfor ;j=0L,size_radi[1]-1  do begin ;radi
endfor ;i=0L,size_thetaR[1]-1  ;thetaR


;draw arrow reference
X0_arrow=+28.
Y0_arrow=-47.

U0=ArrowRef
V0=0.;ArrowRef
mag = SQRT(U0*U0 + V0*V0)

;thetaR(i)= FLOAT(i - L_6lt)*360.0*!DTOR/FLOAT(Xdm[1])
thetaRi=(-45.)*!DTOR
Xmag = V0*cos(thetaRi)-U0*sin(thetaRi)
Ymag = V0*sin(thetaRi)+U0*cos(thetaRi)
costheta = Xmag/mag
if Ymag ge 0 then $
  sintheta = SQRT(1.0 -costheta*costheta) $
else $
  sintheta = - SQRT(1.0 -costheta*costheta)

arwmag = mag * factor
dX_arrow = arwmag * costheta
dY_arrow = arwmag * sintheta

dX_arrow = arwmag * costheta
dY_arrow = arwmag * sintheta

X1_arrow = X0_arrow + dX_arrow
Y1_arrow = Y0_arrow + dY_arrow
;ArrowCol=255. ;white
ARROW, X0_arrow, Y0_arrow, X1_arrow, Y1_arrow  $
, /DATA       $ 
;, /NORMALIZED $
, HSIZE= (!D.X_SIZE / 64.)*0.5 $
;, COLOR=ArrowCol  $ ;index] $
, COLOR=255.  $ ;index] $
;, HTHICK=value] $
;, /SOLID] $
, THICK=value_thick

xyouts,$
; ( (X00+dX-3.2)/X_SIZE ), ( (Y0-0.07)/Y_SIZE) $
 0.70, 0.07 $
, 'V='+STRTRIM( string(ArrowRef,FORMAT='(f5.0)') ,1)+'[m/s]'  $
, charsize=2.0, charthick=2.0, /norm, /noclip

end ;pro draw_arrow_test
