;date 20170804
;purpose: compare with MPA cold plasma density
pro plt_psphere_eqv3_mpa

runYear=13;15
titleYear='20'+STRTRIM( string( (runYear), FORMAT='(i2)'),1 )
print,' titleYear=', titleYear
sw_contourPlot=3L ;1: contour plot; 2: density profiles as function of L; 3: MPA comparison
LValue=[6.19,5.130,3.948,2.925,2.026]
;L Loop
for kL=0,4 do begin 
if kL eq 0 then $ ;geosynchronous
  lpGEOSY=22-1L $       ;      L=6.19546      mlat=66.1310     dL=0.602465
else if kL eq 1 then $
  lpGEOSY=24-1L $       ;      L=5.13083      63.5990     0.462165
else if kL eq 2 then $
  lpGEOSY=27-1L $       ;      L=3.94819      59.5431     0.360240
else if kL eq 3 then $
  lpGEOSY=31-1L $       ;      L=2.92576&nbsp; &nbsp; &nbsp; 53.9256&nbsp; &nbsp; &nbsp;0.226813
else if kL eq 4 then $
  lpGEOSY=37-1L        ;      L=2.02664      44.9684     0.119776
print,kL,' LValue=',LValue[kL],' LpGEOSY=',(LpGEOSY+1)

sw_output2file=0L
sw_output2fileAscii=1L
sw_Plt2Display=0L
sw_dbg=0L
m3ToCm3=1.0e-6

nmp=80L
nlp=93L
NPTS2D=31287L


;read_grid, JMIN_IN,
restore, filename=$
;'~/ipeg/plt/plasma_grid.2xdyn.sav'
;'/Users/naomimaruyama/save/sandbox/ipe/plt.bk20150326/plasma_grid.2xdyn.sav' ;mac
'/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/grid/plt/plasma_grid.2xdyn.sav'

if sw_dbg eq 1 then begin
print, JMIN_IN[0:1]
print, JMAX_IS[0:1]
print, z_km[0]
endif
print,mlat_deg[0:1]

;mlon_deg
mlon_deg=findgen(nmp)*360./FIX(NMP)
if sw_dbg eq 1 then print, mlon_deg

  dt =60L  ;sec
print,'dt=',dt
  pltXsec=900L
  if sw_contourPlot eq 3 then pltXsec = $
;60L ;=ipe_freq_output ;v0 on 201603
900L; v2 on 20170209
print,'pltXsec=',pltXsec

;--
zmax=1.0e+09 ;2.61e+10 ;MAX(z)
zmin=4.79e-9;10 ;MIN(z)
print,'MAX=', zMAX,' MIN=',zMIN
;--
; n read loop
n_read=-1L
nDay=70L
iFileMax=9L
;---
  ut00=$
;518400L ;start_time 00--17ut ;v0
  604800L ;start_time 00--17ut ;v0
;518400+17*3600;583080L  ;17--24ut
  utStart = ut00;+3600*9
  utStop  = utStart+24*3600L*(iFileMax+1) ;17*3600L
  utHrPlt = ut00/3600.
  titleRundir='_ipe_theia_intel_parallel2_93'
;---
if sw_contourPlot eq 3 then begin
;   nMax=1020+1 ;1440L ;v0 
   nMax=97*(iFileMax+1)                   ;v2; on 20170209
   xPlt=fltarr(NMP,nMax)
   yPlt=fltarr(NMP,nMax)
   zPlt=fltarr(NMP,nMax)

if sw_output2fileAscii eq 1 then begin
  lunTmp=2000L
  flnmAsc='equatorialDensity4obana'+'_L'+STRTRIM( string((LValue[kL]), FORMAT='(f4.1)'),1 )+'.YR'+STRTRIM( string( (titleYear), FORMAT='(i4)'),1 )+'.dat'
  print,'opening Ascii file=',flnmAsc
  openw,luntmp,flnmAsc, /GET_LUN

endif
endif
;---iFile loop
for iFile=0, iFileMax do begin
  if sw_dbg eq 1 then print,'iFile=',iFile

;2015
if titleYear eq '2015' then begin
   if iFile eq 0 then $         ;311 
      rundir='1504112979' $            ; 00--24UT 20150311 theia v3 on 20170830
   else if iFile eq 1 then $    ;312 
      rundir='1504133060' $
   else if iFile eq 2 then $    ;313 
      rundir='1504179999' $
   else if iFile eq 3 then $    ;314 
      rundir='1504208116' $
   else if iFile eq 4 then $    ;315 
      rundir='1504263836' $
   else if iFile eq 5 then $    ;316 
      rundir='1504286513' $
   else if iFile eq 6 then $    ;317 
      rundir='1504632658' $
   else if iFile eq 7 then $    ;318 
      rundir='1504681633' $
   else if iFile eq 8 then $    ;319 
      rundir='1504687747' $
   else if iFile eq 9 then $    ;320 
      rundir='1504713868'       ;$
;2013
endif else if titleYear eq '2013' then begin
   if iFile eq 0 then $         ;311 
      rundir='1504129016' $            ; 00--24UT 20130311 theia v3 on 20170830
   else if iFile eq 1 then $    ;312 
      rundir='1504180897' $
   else if iFile eq 2 then $    ;313 
      rundir='1504216211' $
   else if iFile eq 3 then $    ;314 
      rundir='1504265018' $
   else if iFile eq 4 then $    ;315 
      rundir='1504289223' $
   else if iFile eq 5 then $    ;316 
      rundir='1504633830' $
   else if iFile eq 6 then $    ;317 
      rundir='1504684144' $
   else if iFile eq 7 then $    ;318 
      rundir='1504688868' $
   else if iFile eq 8 then $    ;319 
      rundir='1504714210' $
   else if iFile eq 9 then $    ;320 
      rundir='1504804805'
endif ;titleYear eq 
;---
  ut0=0L
  ut0=ut00 - dt ;sec 133.                      ;[ut hr]



  Re_m=6.3712E+06 ;.. Earth radius [meter] 

  x=fltarr(nmp,nlp)
  y=fltarr(nmp,nlp)
  r=fltarr(nmp,nlp)
  z=fltarr(nmp,nlp)
  mlatIN=fltarr(nlp)


;read_plasma, xion
TEST=$
;'tmp20151117'
'mpi20160330v2'
rpath=$
;'~/iper/'+TEST+'/trunk/run/'+rundir+'/'
;'/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/'+TEST+'/run1/'+rundir+'/' ;2015theia
;'/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/'+TEST+'/run/'+rundir+'/' ;2013theia
;'/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/'+TEST+'/run/'+rundir+titleRundir+'/' ;2013theia
'/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/'+TEST+'/run'+STRTRIM( string( (runYear), FORMAT='(i2)'),1 )+'/'+rundir+titleRundir+'/' ;20150311theia
;'/Users/naomimaruyama/save/sandbox/ipe/'+rundir+'/' ;mac

plt_DIR=$
;'~/ipef/'+TEST+'/'+rundir+'/'
;'~/stmp2/'+TEST+'/fig/'+rundir+'/hplus/' ;theia
'/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/fig/plumes/hplus/'
;'/Users/naomimaruyama/sandbox/ipe/fig/'+rundir+'/' ;mac
lun00=0L;h+
lun01=0L;he+
lun02=0L;o+
lun2013=0L
openr,lun00,rpath+'plasma01',/get_lun, /F77_UNFORMATTED ;h+
openr,lun01,rpath+'plasma02',/get_lun, /F77_UNFORMATTED ;he+
openr,lun02,rpath+'plasma00',/get_lun, /F77_UNFORMATTED ;o+
openr,lun2013,rpath+'fort.2013',/get_lun ;, /F77_UNFORMATTED
openr,lun0,rpath+'ut_rec',/get_lun ;, /F77_UNFORMATTED
;read loop
dum=fltarr(NPTS2D,NMP);h+
dum1=fltarr(NPTS2D,NMP);he+
dum2=fltarr(NPTS2D,NMP);o+
ut = 0L
;ut = ut0
record_number=0L



while ( eof(LUN00) eq 0 ) do begin
   readu, lun00,dum ;h+
   readu, lun01,dum1 ;he+
   readu, lun02,dum2 ;o+
   readf, lun0,record_number, ut
;   ut = ut + dt ;[sec]
   n_read = n_read + 1
   print,'rec#',record_number,' ut=', ut,' n_read=',n_read,' iFile=',iFile

   if ut gt utStop then BREAK   ;exit from while loop 

;read sunlon
;   sunlons1 = +0.1030E+01
   readf, LUN2013,sunlons1

   if ut lt utStart then CONTINUE 


; plot only every pltXsec 
   ;if sw_dbg eq 1 then  $
;d print,'ut=',ut,' difut=',(ut-ut00),' MOD difut=',( (ut-ut00)  MOD pltXsec) 
   if (   ((ut-ut00) MOD pltXsec) ne 0 ) then CONTINUE
   ;if sw_dbg eq 1 then  $
   if sw_dbg eq 1 then print,'start plotting ut=',ut


   if ( sw_contourPlot eq 2 AND n_read eq 0 ) then begin
      iwindow=2L
      DEVICE, RETAIN=2, DECOMPOSED=0
      WINDOW,0,XSIZE=700,YSIZE=500
      loadct, n_ldct
   endif

   char_size=1.
   char_thick=1.
   col_max = 255.9999999999999999
   col_min =   0.0000000000000

   axis_color =255.
   n_ldct=0                     ;black+white



   for mp=0,nmp-1 do begin

;calculate MLT, theta
      mltHr    =  mlon_deg[mp]/15.0D0 - sunlons1 * 12.0D0 / !PI  +12.0 ;[hr]
      if ( mltHr lt  0. ) then  mltHr = mltHr MOD 24.
      if ( mltHr ge 24. ) then  mltHr = mltHr - 24.
      mltRad = mltHr*!PI/12.0D0      ;MLT_hr --> THETA[rad]

;note20160519: i am not sure if this is correct?
;shift MLT so that 12MLT on the right!
   ;clockwise 180 deg rotation
      shift_deg= $
 +180.
;- 45. ;original20160519
;      0.
      theta = mltRad - shift_deg/180.*!PI ;(radian)
      for lp=0,nlp-1 do begin
         midpoint = JMIN_IN[lp] + ( JMAX_IS[lp] - JMIN_IN[lp] )/2 -1
         r[mp,lp] = (Z_km[midpoint] * 1.0E+3  + Re_m) / Re_m ;L value
         x[mp,lp] = r[mp,lp] * COS(theta)
         y[mp,lp] = r[mp,lp] * SIN(theta)
         z[mp,lp] = dum[midpoint,mp] ;h+ number density [m-3]     
         mlatIN[lp] = mlat_deg[ JMIN_IN[lp] ]

         if ( sw_contourPlot eq 3 AND lp eq lpGEOSY ) then begin
             if n_read eq 0 AND mp eq 0 then print,' lp=',lp,' mlatIN=',mlatIN[lp]
             if sw_dbg eq 1 then  print,'n_read=',n_read,' mltHr=',mltHr,' mp=',mp

             xPlt[mp,n_read]=mltHr
             yPlt[mp,n_read]=FIX(nDay) + (ut-ut00)/86400.
             zPlt[mp,n_read]=dum[midpoint,mp]*m3ToCm3 ;h+ number density [cm-3]  


if (sw_output2fileAscii eq 1) AND (yPlt[mp,n_read] ge 73.00) then $
  printf,lunTmp,yPlt[mp,n_read],xPlt[mp,n_read], (dum[midpoint,mp]*m3ToCm3) $ ;h+
,(dum1[midpoint,mp]*m3ToCm3) $ ;he+
,(dum2[midpoint,mp]*m3ToCm3) $ ;o+
, FORMAT='(F7.3,F7.2,3e12.4)'

        
         endif                    ;sw_contourPlot=3   

      endfor                    ; lp=0,nlp-1 do begin


      

      
      if (sw_contourPlot eq 2) then begin

         m2cm=1.0e+02

         xmax=5.
         xmin=2.
         if (record_number eq 1 AND mp eq 1 ) then $
            plot,r[mp,*],(z[mp,*]*1.0E-6) $
;         ,xrange=[2.    ,5.    ], xstyle=1 $
                 ,xrange=[xmin    ,xmax    ], xstyle=1 $
                 ,yrange=[1.e+00,1.e+06], ystyle=1 $
                 ,/YLOG $
                 ,title='h+ number density'  $
                 ,linestyle = 0 $
                 ,color=axis_color $
                 ,charsize=char_size $
                 ,/NODATA


;debug print

;      if ( sw_debug eq 1 ) then begin
         if ( mp eq 0 ) then begin
            for lp=1,nlp-1 do begin
               if lp eq 1 then print, 'lp,  r,  mlatIN, dlat, dr'
               print, lp,r[mp,lp],mlatIN[lp],(mlatIN[lp-1]-mlatIN[lp]), (r[mp,lp-1]-r[mp,lp])
                                ; if ( r[mp,lp] ge 2. AND r[mp,lp] lt 5 ) then begin
                                ;    print, 'check r&z', r[mp,lp],(z[mp,lp]*1.0E-6),mlatIN[lp]
                                ; endif
            endfor              ;lp
         endif                  ;( mp eq 0 ) then begin
;      endif
;      STOP

;      
         n_ldct=39              ;black+white
         loadct, n_ldct
   
         oplot,r[mp,0:nlp-1],(z[mp,0:nlp-1]*1.0E-6) $
               ,linestyle=0 $  
               ,color = mlt * col_max / (2.*!PI)


 



      
      
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
;choose color
 text_color=col_min 

if ( sw_contourPlot eq 1 ) then begin 
   iwindow=1L
   DEVICE, RETAIN=2, DECOMPOSED=0
   WINDOW,0,XSIZE=500,YSIZE=500
   loadct, n_ldct

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
            , LINESTYLE = 0
  loadct, n_ldct


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
print,'ut before xyout', ut,  ut/3600.,  (ut/3600.-utHrPlt)
xyouts, 0.015, 0.02  $
,'UT [hrs]='+STRTRIM( string( (ut/3600.-utHrPlt), FORMAT='(F7.3)'),1 ) $
, charsize =1.5, charthick=1.5 $
, /norm, /noclip


xyouts, 0.80, 0.02  $
,rundir $
, charsize =0.9, charthick=0.9 $
        , /norm, /noclip

    if ( sw_output2file eq 1 ) then begin
       filename_png=plt_DIR+'hp_mgeq_ut'+STRTRIM( string((ut/3600.-utHrPlt), FORMAT='(F6.2)'),1 )+rundir+'v1_180.png'
       if sw_dbg eq 1 then print, filename_png
       output_png, filename_png
    endif                       ; ( sw_output2file eq 1 ) then begin

 endif  else if (sw_contourPlot eq 3 AND n_read eq nMax-1) then begin 

;  MPA comparison    
if sw_Plt2Display eq 1 then begin
    print,"plotting MPA comparison", ut, ut00, (ut-ut00),nMax

    iwindow=1L
    DEVICE, RETAIN=2, DECOMPOSED=0
    WINDOW,0,XSIZE=1000,YSIZE=1000
    n_ldct=39
    loadct, n_ldct

    n_levels=100
    zmax=100.
    zmin=0.
    x_max=24.
    x_min=0.
    y_max=80.
    y_min=70.
    
    
    char_size=1.
    char_thick=1.
    col_max = 255.9999999999999999
    col_min =   0.0000000000000
    axis_color =255.
;    n_ldct=0                    ;black+white
    text_color=col_min;col_min 
    
    print,'MAX=', MAX(zPlt),MIN(zPlt)

    loadct,39
    contour,zplt,xplt,yplt $
           ,/irregular $
           ,/fill $
           , levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
           , xrange=[X_min,X_max], /xstyle  $
           , yrange=[Y_min,Y_max], /ystyle  $
           ,XTITLE = 'MLT', YTITLE = 'DOY' $ 
           ,TITLE = 'LANL/MPA-IPE Comparison: '+titleYear $
;, POSITION=[X0 , Y0 , X1 , Y1 ] $
           , COLOR=255 $;text_color $
           , charsize = char_size, charthick = char_thick ;$
;, MAX_VALUE= MAX_xymin ;$
;           ,/NODATA
;---



title=' ';colorbar title'
font=1   ;True-Type: 1.
; X_SIZE=29.7
; Y_SIZE=20.5 
; X0 = 2.2
; Y0 = 2.6
; Y_cmargin=1.0
;dYc = 0.28
charsize_colorbar=2.
format_colorbar='(F7.0)' ;(E9.1)'
;position=[X0/X_SIZE, (Y1+Y_cmargin)/Y_SIZE, X1/X_SIZE, (Y1+Y_cmargin+dYc)/Y_SIZE]
position=[0.095, 0.97, 0.92, 0.98]
COLORBAR, BOTTOM=bottom, CHARSIZE=charsize_colorbar, COLOR=color, DIVISIONS=divisions $
        , FORMAT=format_colorbar, POSITION=position, MAXRANGE=zmax, MINRANGE=zmin $
        , NCOLORS=ncolors, TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
          , _EXTRA=extra, INVERTCOLORS=invertcolors,  TICKNAMES=ticknames

    if ( sw_output2file eq 1 ) then begin
       filename_png=plt_DIR+'hp_MPA_lp'+STRTRIM( string((lpGEOSY+1), FORMAT='(i3)'),1 )+'_'+rundir+'.png'
       if sw_dbg eq 1 then print, filename_png
       output_png, filename_png
    endif                       ; ( sw_output2file eq 1 ) then begin

endif else if sw_Plt2Display eq 0 then print,'plotting contour skipped'

 endif                          ;sw_contourPlot   


endwhile                        ;( eof(LUN00) eq 0 ) do begin

free_Lun,lun00
free_Lun,lun01
free_Lun,lun02
free_Lun,lun0
free_Lun,lun2013
endfor ;iFile=0, 1 do begin
endfor ;kL=0,4 do begin 
print, 'plt_psphere_eqv3_MPA finished!'
end                             ;pro plt_psphere_eqv3_MPA
