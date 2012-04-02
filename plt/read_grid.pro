pro read_grid,LUN,JMIN_IN,JMAX_IS,Z_km,mlat_deg,sw_debug

     mp=0L

if sw_debug eq 1 then  print,'size jmin',SIZE(jmin_in)
if sw_debug eq 1 then  print,'size jmax',SIZE(jmax_is)




if sw_debug eq 1 then  print,'size z_km',SIZE(Z_km)
if sw_debug eq 1 then  print,'size mlat_deg ',SIZE(mlat_deg)

NMP_all=0L
NLP_all=0L
NPTS2D_dum=0L
; read plasma_grid
    readu, LUN[1], NMP_all
print,'NMP_all',NMP_all
    readu, LUN[1], NLP_all
print,'NLP_all',NLP_all
    readu, LUN[1], NPTS2D_dum
print,'NPTS2D_dum',NPTS2D_dum


lpj=0L

    readu, LUN[1], JMIN_IN ;(1:NLP_all)
print,lpj,'JMIN_IN(1)',JMIN_IN[lpj]
    readu, LUN[1], JMAX_IS ;(1:NLP_all)
print,'JMAX_IS(1)',JMAX_IS[lpj]

lpj=129L
ipts=JMIN_IN(lpj)-1

    dum = fltarr( NMP_all+1 ) ;!rad
    readu, LUN[1], dum ;( 1: NMP_all+1 ) !rad
print,'mlon[deg]',dum[mp]*180./!PI

    dum = Z_km
    readu, LUN[1], dum ;(     1:NPTS2D_dum) !meter
     Z_km     = dum * 1.0E-3
print,'z_km',z_km[ipts]
    dum  = mlat_deg
    readu, LUN[1], dum ;(    1:NPTS2D_dum) !rad
     mlat_deg = ( !PI*0.50 - dum ) * 180.0 / !PI
print,'mlat_deg',mlat_deg[ipts]

    dum=fltarr(NPTS2D_dum, NMP_all)
    readu, LUN[1], dum
print,'GCOLAT-deg',90.-dum[ipts, mp]*180./!PI

    dum=fltarr(NPTS2D_dum, NMP_all)
    readu, LUN[1], dum
print,'GLON_rad-deg',dum[ipts, mp]*180./!PI  

;     readu, LUN[1], JMIN_IN,JMAX_IS,Z_meter,GL_rad


if sw_debug eq 1 then  print, 'plasma0: IN=',JMIN_IN[0],' IS=',JMAX_IS[0],Z_km[0],mlat_deg[0]

;STOP ;debug
END ;PRO read_grid
