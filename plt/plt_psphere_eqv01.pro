;underconstruction
;comlat[deg] vs. R???
;v01: add EXB drift
;i should also add polar cap potential variations
pro plt_psphere_eqv01

titleYear='2013'
sw_LNewZCircle=0;1
sw_SedPlumeFluxTube=0;1

;0:no arrow
;1:arrow read in vexbup
;2:calculate Vr
sw_arrow=2L
;###CHANGE
titlePlot= $
;'2013-3-17 quiet';run2
'2013-3-17 storm' ;run'
;'2015-3-17 storm' ;run1
sw_contourPlot=1L ;1: contour plot; 2: density profiles as function of L
sw_output2file=1L
sw_dbg=0L
  dt =900L  ;sec
  pltXsec=900L
;###CHANGE
  ut00=$
;518400L ;start_time 00--17ut
;579600;518400+17*3600; ;17--24ut
583080
;604800L ;start_time 00--24ut ;v2
  utStart=ut00
  utStop=utStart;+24*3600;
  utHrPlt=ut00/3600.
titleRundir='_ipe_theia_intel_parallel2_93'
;###CHANGE
rundir=$
;'ipe_S_32328' ;original
;'ipe_S_25827' ;depeleted
;'ipe_S_26060'; ;transport only
;'1461312397_ipe_theia_intel_parallel2_93';20130317 00--17 before dep
'1461343191' ;20130317 17--24 ;mac
;'1461417243_ipe_theia_intel_parallel2_93';.00_17UT20150317';mac
;'1461436195_ipe_theia_intel_parallel2_93'; 17ut
;'1463696937_ipe_theia_intel_parallel2_93' ;2013 theia after dep
;'1462618349_ipe_theia_intel_parallel2_93' ;2013 00--17ut quiet
;'1462657622_ipe_theia_intel_parallel2_93' ;17ut quiet
;'1486587459'; 00--24UT 20130317 theia v2 on 20170209
;---
  ut0=0L
  ut0=ut00 - dt ;sec 133.                      ;[ut hr]

;--
;zmax=1.0e+09 ;2.61e+10 ;MAX(z)
;zmin=4.79e-9;10 ;MIN(z)
zmax=100. ;for L>=4
zmin=0.
print,'MAX=', zMAX,zMIN
;--
  nmp=80L
  nlp=93L
  NPTS2D=31287L
  Re_m=6.3712E+06 ;.. Earth radius [meter] 

  x=fltarr(nmp,nlp)
  y=fltarr(nmp,nlp)
  r=fltarr(nlp)
  z=fltarr(nmp,nlp)

;read_grid, JMIN_IN,
restore, filename=$
'/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/grid/plt/plasma_grid.2xdyn.sav';theia
;'/Users/naomimaruyama/sandbox/ipe/plt.bk20150326/plasma_grid.2xdyn.sav';mac

if sw_dbg eq 1 then begin
print, JMIN_IN[0:1]
print, JMAX_IS[0:1]
print, z_km[0]
endif

;mlon_deg
mlon_deg=findgen(nmp)*360./FIX(NMP)
if sw_dbg eq 1 then print, mlon_deg

;read_plasma, xion
TEST=$
;'tmp20151117'
'mpi20160330v2'
;'depletedFlux20160512'
;###CHANGE
hpath='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/'
rpath=$
;'~/iper/'+TEST+'/trunk/run/'+rundir+'/'
;  hpath+TEST+'/run/'+rundir+'/' ;theia 2013 storm
;  '~/stmp2/'+TEST+'/run1/'+rundir+'/' ;theia 2015
;  '~/stmp2/'+TEST+'/run2/'+rundir+'/' ;theia quiet
;'/Users/naomimaruyama/sandbox/ipe/'+rundir+'/' ;mac
'/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/'+TEST+'/run/'+rundir+titleRundir+'/' ;2013theia
plt_DIR=$
;'~/ipef/'+TEST+'/'+rundir+'/'
;   '~/stmp2/'+TEST+'/fig/'+rundir+'/hplus/' ;theia
;'/Users/naomimaruyama/sandbox/ipe/fig/'+rundir+'/' ;mac
'/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/fig/plumes/hplus/'
lun00=0L
lun17=0L
lun16=0L
lun18=0L
lun2013=0L
openr,lun00,rpath+'plasma01',/get_lun, /F77_UNFORMATTED ;h+
openr,lun17,rpath+'plasma17',/get_lun, /F77_UNFORMATTED ;VEXBe
openr,lun16,rpath+'plasma16',/get_lun, /F77_UNFORMATTED ;VEXBup
openr,lun18,rpath+'plasma18',/get_lun, /F77_UNFORMATTED ;VEXBth
openr,lun2013,rpath+'fort.2013',/get_lun ;, /F77_UNFORMATTED
openr,lun0,rpath+'ut_rec',/get_lun ;, /F77_UNFORMATTED
;read loop
dum=fltarr(NPTS2D,NMP)
vexbe=fltarr(NLP,NMP);plasma17
vexbup=fltarr(NLP,NMP);plasma16
vexbth=fltarr(NLP,NMP);plasma17
ut = 0L
;ut = ut0
min_record_number=1079L
record_number=0L
while ( eof(LUN00) eq 0 ) do begin
   readu, lun00,dum
   readu, lun17,vexbe
   readu, lun16,vexbup
   readu, lun18,vexbth
   readf, lun0,record_number, ut
   print,'rec#',record_number,' ut=', ut
;   ut = ut + dt ;[sec]
   if ut gt utStop then STOP 

;read sunlon
;   sunlons1 = +0.1030E+01
   readf, LUN2013,sunlons1

   if ut lt utStart then CONTINUE 


; plot only every pltXsec 
   print,ut,'difut=', (ut-ut00), pltXsec, ( (ut-ut00)  MOD pltXsec) 
   if (   ((ut-ut00) MOD pltXsec) ne 0 ) then CONTINUE
   print,'start plotting ut=',ut


   if ( sw_contourPlot eq 2 AND record_number eq min_record_number ) then begin
      iwindow=2L
      DEVICE, RETAIN=2, DECOMPOSED=0
      WINDOW,0,XSIZE=700,YSIZE=500
   endif

   char_size=1.
   char_thick=1.
   col_max = 255.9999999999999999
   col_min =   0.0000000000000

   axis_color =255.
   n_ldct=0                     ;black+white
   loadct, n_ldct

   mlt=fltarr(nmp)
   for mp=0,nmp-1 do begin

;calculate MLT, theta
      mlt[mp]    =  mlon_deg[mp]/15.0D0 - sunlons1 * 12.0D0 / !PI  +12.0 ;[hr]

      ;if mp eq 56 then print, 'mlt', mlt[mp]-24

      if ( mlt[mp] lt  0. ) then  mlt[mp] = mlt[mp] MOD 24.
      if ( mlt[mp] ge 24. ) then  mlt[mp] = mlt[mp] - 24.
      mlt[mp] = mlt[mp]*!PI/12.0D0      ;MLT_hr --> THETA[rad]

;note20160519: i am not sure if this is correct?
;shift MLT so that 12MLT on the right!
   ;clockwise 180 deg rotation
      shift_deg= $
+180. ;shift 180deg clockwise
;- 45. ;original20160519
;      0.  ;tried 20160519
      mlt[mp] = mlt[mp] + shift_deg/180.*!PI ;(radian)
      for lp=0,nlp-1 do begin
         midpoint = JMIN_IN[lp] + ( JMAX_IS[lp] - JMIN_IN[lp] )/2 -1
         r[lp] = (Z_km[midpoint] * 1.0E+3  + Re_m) / Re_m ;L value
         x[mp,lp] = r[lp] * COS(mlt[mp])
         y[mp,lp] = r[lp] * SIN(mlt[mp])
         z[mp,lp] = dum[midpoint,mp]*1.E-6 ;h+ number density [cm-3]     
      endfor                         ; lp=0,nlp-1 do begin




   if (sw_contourPlot eq 2) then begin

      m2cm=1.0e+02

      xmax=5.
      xmin=2.
      if (record_number eq min_record_number AND mp eq 1 ) then $
    plot,r[*],(z[mp,*]*1.0E-6) $
;         ,xrange=[2.    ,5.    ], xstyle=1 $
              ,xrange=[xmin    ,xmax    ], xstyle=1 $
     ,yrange=[1.e+00,1.e+06], ystyle=1 $
     ,/YLOG $
  ,title=titlePlot+' h+ number density'  $
  ,linestyle = 0 $
  ,color=axis_color $
  ,charsize=char_size $
  ,/NODATA


;debug print

      if ( sw_debug eq 1 ) then begin
         for lp=0,nlp-1 do begin
            print, mp,lp,r[lp]
            if ( r[lp] ge 2. AND r[lp] lt 5 ) then begin
               print, 'check r&z', r[lp],(z[mp,lp]*1.0E-6)
            endif
         endfor                 ;lp
      endif
      
      n_ldct=39                  ;black+white
      loadct, n_ldct
   
      oplot,r[0:nlp-1],(z[mp,0:nlp-1]*1.0E-6) $
            ,linestyle=0 $  
            ,color = mlt[mp] * col_max / (2.*!PI)


 endif                          ;(sw_contourPlot eq 2) then begin
    
endfor                          ;mp=0,nmp-1 do begin    




n_levels=100
X_max=+7.0
X_min=-X_max
Y_max=X_max
Y_min=-Y_max



;MAX_xymin =

;n_ldct=8 ;green
n_ldct=1 ;blue
loadct, n_ldct
;choose color
 text_color=col_min 

if ( sw_contourPlot eq 1 ) then begin 
iwindow=1L
DEVICE, RETAIN=2, DECOMPOSED=0
WINDOW,0,XSIZE=500,YSIZE=500

print,'MAX=', MAX(z),MIN(z)

contour,z,x,y $
,/irregular $
,/fill $
, levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
, xrange=[X_min,X_max], /xstyle  $
, yrange=[Y_min,Y_max], /ystyle  $
,XTITLE = 'X', YTITLE = 'Y' $ 
,TITLE = ' ' $
;, POSITION=[X0 , Y0 , X1 , Y1 ] $
, COLOR=text_color $
, charsize = char_size, charthick = char_thick  ;$
;, MAX_VALUE= MAX_xymin ;$
;,/NODATA


;add arrow
if sw_arrow ge 1 then begin
   u=fltarr(nmp,nlp)
   v=fltarr(nmp,nlp)
;assign u,v
   for lp=0,nlp-1 do begin
      for mp=0,nmp-1 do begin


         u[mp,lp] = vexbe[lp,mp]*1.0e-3
         v[mp,lp] = vexbup[lp,mp]*1.0e-3
         if sw_arrow eq 2 then begin
            calculate_vr, vr,ve,mp,lp,VEXBTH,VEXBE,z_km,mlat_deg,jmin_in,jmax_is
;print,mp,lp,'Vr=',vr,' VEXBup=',vexbup[lp,mp]
if lp eq 21 then print,mp,lp,'Vr[m/s]=',vr,' VEXBup',vexbup[lp,mp],' Ve=',ve,' VEXBe=',vexbe[lp,mp]
            u[mp,lp]=ve*1.0e-3 ;m/s==>km/s
            ;v[mp,lp]=vr*1.0e-3 ;m/s==>km/s
         endif ;sw_arrow
      endfor                    ;mp
   endfor                       ;lp
;plot VEXB
;mlt[rad]
;comlat[deg] measured from North pole-->R
;rim_lat
; draw_arrow_test, u, v, mlt, comlat, rim_lat, sw_debug $
;arrow ref location
   X0_arrow=(6.6-6.6) * 0.563       ;+28.
   Y0_arrow=(6.6-6.6) * 0.9451 * (-1.) ;-47.
   draw_arrow_test3, u, v, mlt, r,  sw_debug $
                    ,x0_arrow,y0_arrow
endif                           ;_sw_arrow eq 1 then begin


; plot earth
;---earth start
; potential inside the earth should not be plotted!!!
; create the earth artificially
; zz=dblarr(isize,jsize)
jsize_e = 200L
isize_e = 200L
xmin_e=fltarr(isize_e,jsize_e)
ymin_e=fltarr(isize_e,jsize_e)
zz_e  =fltarr(isize_e,jsize_e)
rxy_max = 1.0
xmin_e0 = -1.0
ymin_e0 = -1.0

for j = 1-1, jsize_e-1  do begin
  for i = 1-1, isize_e-1  do begin

    xmin_e[i,j] =  xmin_e0 + i*rxy_max/isize_e*2.
    ymin_e[i,j] =  ymin_e0 + j*rxy_max/jsize_e*2.
    rxy = SQRT( xmin_e[i,j]*xmin_e[i,j] + ymin_e[i,j]*ymin_e[i,j] )    

    if ( rxy le rxy_max ) then begin
       if ( xmin_e[i,j] ge 0. ) then  zz_e[ i,j] = +rxy_max*0.7 $
       else                           zz_e[ i,j] = -rxy_max*0.7

    endif else begin ;if ( rxy > rxy_max ) then begin
       xmin_e[i,j]=0.0
       ymin_e[i,j]=0.0
       zz_e[  i,j]=0.0
    endelse

  endfor
endfor
zmax_e =  rxy_max
zmin_e = -rxy_max
;---earth end

loadct, 0 ;blue/white linear
contour,zz_e,xmin_e,ymin_e $
,/irregular $
,/fill $
, levels=findgen(n_levels)*(zmax_e-zmin_e)/float(n_levels-1) +zmin_e $
;, xrange=[X_min,X_max], /xstyle  $
;, yrange=[Y_min,Y_max], /ystyle  $
;, XTITLE = 'XMIN' , YTITLE = Y_TITLE $
;, TITLE = ' ' $
;, POSITION=[ X0/X_SIZE, Y0/Y_SIZE, X1/X_SIZE, Y1/Y_SIZE] $
;tmp020810: , COLOR=text_color $
;, charsize = char_size, charthick = char_thick  $
, /NOERASE $
, /overplot

; add circles along L values
  loadct, 0
  col_cir = col_max  ;white
  ;col_cir = col_min ;black

  theta_rad = findgen(360)*1.0D0 /180. * !PI ;[deg]-->[rad]
  for k=1,5,2  do begin

;if k ne 3 OR k ne 5 then CONTINUE 

     radius_circle = findgen(360)*0.0 + 1.0*FLOAT(k+1)

;back ground black circle
    oplot , /polar , radius_circle , theta_rad  $
            , THICK=0.8 $
            , COLOR =0. $;black  col_cir $
            , LINESTYLE = 0
;over plotting white circle
    oplot , /polar , radius_circle , theta_rad  $
            , THICK=0.8 $
            , COLOR =255. $     ;white  col_cir $
            , LINESTYLE = 2
  endfor
;geo synchronous
  radius_circle = findgen(360)*0.0 + 6.6
  oplot , /polar , radius_circle , theta_rad  $
            , THICK=1.0 $
            , COLOR =255. $     ;white  col_cir $
;            , LINESTYLE = 0
            , LINESTYLE = 1
  loadct, n_ldct




;circle NewZ land flux tube: lpSed=31L--> L=2.74343
if ( sw_LNewZCircle eq 1 ) then begin
  n_ldct_sv = n_ldct
  loadct, 39

  Lx=2.74343
  radius_circle = findgen(360)*0.0 + Lx
    oplot , /polar , radius_circle , theta_rad  $
            , THICK=0.6 $
            , COLOR =250. $     ;red
            , LINESTYLE = 2
endif ;( sw_LNewZCircle eq 1 ) then begin

;following a specific flux tube:
if ( sw_SedPlumeFluxTube eq 1 ) then begin
  mpSed=3L
  lpSed=31L
  oplot,   x[mpSed,lpSed-1:lpSed], y[mpSed,lpSed-1:lpSed]  $
            , THICK=9.0, LINESTYLE = 0  $
            , COLOR = 250. ;red

  oplot,   x[mpSed:mpSed+1,lpSed], y[mpSed:mpSed+1,lpSed]  $
            , THICK=9.0, LINESTYLE = 0  $
            , COLOR = 250. ;green
  loadct, n_ldct_sv
endif ;( sw_PlumeFluxTube eq 1 ) then begin
;

;colorbar
title=' ';colorbar title'
font=1   ;True-Type: 1.
; X_SIZE=29.7
; Y_SIZE=20.5 
; X0 = 2.2
; Y0 = 2.6
; Y_cmargin=1.0
;dYc = 0.28
charsize_colorbar=2.
format_colorbar='(E9.1)'
;position=[X0/X_SIZE, (Y1+Y_cmargin)/Y_SIZE, X1/X_SIZE, (Y1+Y_cmargin+dYc)/Y_SIZE]
position=[0.08, 0.96, 0.92, 0.98]
COLORBAR, BOTTOM=bottom, CHARSIZE=charsize_colorbar, COLOR=color, DIVISIONS=divisions $
        , FORMAT=format_colorbar, POSITION=position, MAXRANGE=zmax, MINRANGE=zmin $
        , NCOLORS=ncolors, TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
          , _EXTRA=extra, INVERTCOLORS=invertcolors,  TICKNAMES=ticknames


loadct,0
if sw_dbg eq 1 then  print,'ut before xyout', ut,  ut/3600.,  (ut/3600.-utHrPlt)
xyouts, 0.015, 0.02  $
,'UT [hrs]='+STRTRIM( string( (ut/3600.-utHrPlt), FORMAT='(F7.3)'),1 ) $
, charsize =1.5, charthick=1.5 $
, /norm, /noclip


xyouts, 0.8, 0.026  $
,titlePlot $
, charsize =1., charthick=1. $
        , /norm, /noclip
xyouts, 0.68, 0.008  $
,rundir $
, charsize =0.8, charthick=0.8 $
        , /norm, /noclip


utDisp=ut/3600.-utHrPlt
if utDisp lt 10. then $
   stringUt='0'+STRTRIM( string(utDisp, FORMAT='(F6.2)'),1 ) $
else $; if utDisp ge 10.
   stringUt=STRTRIM( string(utDisp, FORMAT='(F6.2)'),1 )

if ( sw_output2file eq 1 ) then begin
   filename_png= $ 
;if mpSed is set
;plt_DIR+'hp_mgeq_ut'+stringUt+rundir+'mp'+STRTRIM( string(mpSed, FORMAT='(i2)'),1 )+'.png'
plt_DIR+'hp_mgeq_ut'+stringUt+rundir+'.png'
if sw_dbg eq 1 then print, filename_png
   output_png, filename_png
endif; ( sw_output2file eq 1 ) then begin

endif ;( sw_contourPlot eq 1 ) then begin

endwhile                        ;( eof(LUN00) eq 0 ) do begin

free_Lun,lun00
free_Lun,lun0
free_Lun,lun2013
free_Lun,lun17
free_Lun,lun16
print, 'plt_psphere_eqv01 finished!'
end                             ;pro plt_psphere_eqv0
