pro ctr_lon_lat $
,JMIN_IN,JMAX_IS,Z_km,mlat_deg $
;,je_3d $
, XIONN_m3, TE_TI_k $
, XIONV_ms1 $
,UT_hr, plot_DIR $
,n_read $
,sw_output2file $
,glon_deg,glat_deg,sw_frame $ ;=0 ;mag; 1geo
,fac_window ,TEST $
,sw_debug


for VarType=0,0,5 do begin

print,'VarType=',VarType
ht_plot =350.;95.;130.;250.;91.;380.;180.;900.;380.;280.;380.00 ;[km]
;ht_plot =500.;280.;380.00 ;[km]
;ht_plot =600.00 ;[km]


sw_range=1L
nano=1.0E-9
unit=['[m-3]','K','K','[m-3]','[m-3]','[m-3]','m/s'];[nA/m2]'];[nA/m2]'

 sw_plot_grid=0 ;1:retangular, 2:polar
; get n_read_max, ny_max
size_result = SIZE(JMIN_IN)
if ( sw_debug eq 1 ) and ( n_read eq 0 ) then  print,'NLP',size_result
NLP = size_result[1]
ny_max =NLP*2
if  ( sw_debug eq 1 ) and ( n_read eq 0 ) then $
   print,'NLP',NLP,'ny_max',ny_max

size_result = SIZE(XIONN_m3) ;je_3d)
if  ( sw_debug eq 1 ) and ( n_read eq 0 ) then  print,'NMP',size_result
ISPEC=size_result[1]
NPTS2D=size_result[2]
NMP = size_result[3]
nx_max =size_result[3]    ;NMP=80
if  ( sw_debug eq 1 ) and ( n_read eq 0 ) then $
  print,'NMP=',NMP,' nx_max=',nx_max,' ISPEC=',ispec,' NPTS2D=',NPTS2D


mlon_deg_min=0. ;default
dmlon=4.50 ;deg
mlon_deg=findgen(NMP)*dmlon
for mp=0,nmp-1 do begin
if ( mlon_deg[mp] lt mlon_deg_min ) then  mlon_deg[mp]=mlon_deg[mp]+360.
endfor

plot_zz=fltarr(nx_max,ny_max)
plot_yy=fltarr(nx_max,ny_max)
plot_xx=fltarr(nx_max,ny_max)

VarTitle=['Ne','Te','Ti','O+','NO+','O2+','Vpar_o+']
;VarTitle=['Ne','Te','Ti','N(NO+)','Vpar_o+']
;VarTitle=['Ne','Te','Ti','N(O2+)','Vpar_o+']
;VarTitle=['Ne','Te','Ti','N(N2+)','Vpar_o+']

for mp=0,NMP-1 do begin
for lp=0,NLP-1 do begin

  in = JMIN_IN[lp]-1L
  is = JMAX_IS[lp]-1L
  midpoint=IN + ( IS - IN )/2
;NH  
  istep=+1
  for i=in,midpoint, istep  do begin
    if ( z_km[i] le ht_plot ) AND ( z_km[i+istep] gt ht_plot ) then begin
;print,'NH: mp',mp,' lp',lp,' i',i,' z_km',z_km[i] , ht_plot,XIONN_m3[VarType,i-20:i+20,mp] 
      if ( VarType eq 0 ) then $ 
        for jth=0,ISPEC-1 do  plot_zz[mp,lp]=plot_zz[mp,lp]+XIONN_m3[jth,i,mp]*1.0E-10 $
      else if ( VarType eq 1 ) then $
        plot_zz[mp,lp] = TE_TI_k[3-1,i,mp] $;[K]
      else if ( VarType eq 2 ) then $
        plot_zz[mp,lp] = TE_TI_k[1-1,i,mp] $;[K]
      else if ( VarType eq 3 ) then $
        plot_zz[mp,lp] = XIONN_m3[1-1,i,mp]*1.0E-10  $
      else if ( VarType eq 4 ) then $
        plot_zz[mp,lp] = XIONN_m3[4,i,mp]*1.0E-10  $ ;no+
      else if ( VarType eq 5 ) then $
        plot_zz[mp,lp] = XIONN_m3[5,i,mp]*1.0E-10  $ ;o2+
      else if ( VarType eq 6 ) then $
        plot_zz[mp,lp] = XIONV_ms1[1-1,i,mp] ;V//o+[m/s]
 
       if ( sw_frame eq 0 ) then begin ;magnetic
          plot_yy[mp,lp] = mlat_deg[i]
          plot_xx[mp,lp] = mlon_deg[mp]
;          if (  mlon_deg[mp] ge 180.) then plot_xx[mp,lp]=mlon_deg[mp]-360.
       endif else if ( sw_frame eq 1 ) then begin ;geographic
          plot_yy[mp,lp] = glat_deg[i,mp]
          plot_xx[mp,lp] = glon_deg[i,mp]
;          if (  glon_deg[i,mp] ge 180.) then plot_xx[mp,lp]=glon_deg[i,mp]-360.
       endif ;( sw_frame eq 0 ) then begin ;magnetic
       BREAK 
    endif
  endfor ;i

;SH  
  lps = NLP-1 + (NLP-1-lp)
  istep=-1 
  for i=is,midpoint, istep  do begin
    if ( z_km[i] le ht_plot ) AND ( z_km[i+istep] gt ht_plot )  then begin

;print,'SH: mp',mp,' lp',lp,' i',i,' z_km',z_km[i] , ht_plot,XIONN_m3[VarType,i-20:i+20,mp] 

;       plot_zz[mp,lps] = XIONN_m3[VarType,i,mp] ;[m-3]  je_3d[VarType,i,mp]/nano
      if ( VarType eq 0 ) then $ 
        for jth=0,ISPEC-1 do  plot_zz[mp,lps]=plot_zz[mp,lps]+XIONN_m3[jth,i,mp]*1.0E-10 $
      else if ( VarType eq 1 ) then $
        plot_zz[mp,lps] = TE_TI_k[3-1,i,mp] $;[K]
      else if ( VarType eq 2 ) then $
        plot_zz[mp,lps] = TE_TI_k[1-1,i,mp] $;[K]
      else if ( VarType eq 3 ) then $
        plot_zz[mp,lps] = XIONN_m3[1-1,i,mp]*1.0E-10 $
      else if ( VarType eq 4 ) then $
        plot_zz[mp,lps] = XIONN_m3[4,i,mp]*1.0E-10 $
      else if ( VarType eq 5 ) then $
        plot_zz[mp,lps] = XIONN_m3[5,i,mp]*1.0E-10 $
      else if ( VarType eq 6 ) then $
        plot_zz[mp,lps] = XIONV_ms1[1-1,i,mp] ;V//o+[m/s]

       if ( sw_frame eq 0 ) then begin ;magnetic
          plot_yy[mp,lps] = mlat_deg[i]
          plot_xx[mp,lps] = mlon_deg[mp]
;          if (  mlon_deg[mp] ge 180.) then plot_xx[mp,lps]=mlon_deg[mp]-360.
       endif else if ( sw_frame eq 1 ) then begin ;geographic
          plot_yy[mp,lps] = glat_deg[i,mp]
          plot_xx[mp,lps] = glon_deg[i,mp]
;          if (  glon_deg[i,mp] ge 180.) then plot_xx[mp,lps]=glon_deg[i,mp]-360.
       endif ;( sw_frame eq 0 ) then begin ;magnetic
       BREAK
    endif
 endfor                         ;i

endfor ;lp=0,NLP-1 do begin
endfor ;mp=0,NMP-1 do begin


;(2) when time = time_max 
; plotting
if ( sw_range eq 1 ) then begin
      if ( VarType eq 0 ) then begin 
;f107-180
        zmin=1.
;        zmin=1.e+10
;        zmax=7.5e+11
;        zmax=2.5e+1 ;E-region
        zmax=2.5e+2 ;F-region
;        zmax=5.58e+12
;92km
;        zmax=2.00e+11
;        zmin=1.00e+04
;f107-80
;        zmin=1.e+10
;        zmax=1.77e+12
      endif else if ( VarType ge 1 ) AND (VarType le 2 ) then begin 
        zmin=440.
        zmax=640.
      endif else if ( VarType ge 3 ) AND ( VarType le 5 ) then begin 
;        zmin=0.
;        zmax=8.
        zmin=1.
        zmax=2.5e+1
      endif else if ( VarType eq 6 ) then begin 
        zmin=-200.
        zmax=+200.
      endif

	where_result=where ( plot_zz gt zmax, count ) 
;	print, '1. where_result', where_result,'count',count
	for i=0,count-1 do $
	  plot_zz[ where_result[i] ]=zmax
;  	where_result=where ( plot_zz gt zmax, count ) 
;	print, '2. where_result', where_result,'count',count

endif else if ( sw_range eq 0 ) then begin
zmax = max(plot_zz)
zmin = min(plot_zz)
endif
;if ( sw_debug eq 1 ) then   
print,'max=',MAX(plot_zz),' min=',MIN(plot_zz)

MAX_xymin=zmax ;1.0E+3
n_levels=100L

X_max=+360.;mlon_deg_min+360.
X_min=  00.;mlon_deg_min ;  0.0
Y_max=+90.;60.;90.;60.0
Y_min=-Y_max
if ( sw_frame eq 0 ) then $
  MAG_GEO='magnetic' $
else if ( sw_frame eq 1 ) then $
  MAG_GEO='geographic'

X_TITLE=MAG_GEO+' longitude[deg]'
Y_TITLE=MAG_GEO+' latitude[deg]'

text_color=255.
char_size=1.0
char_thick=1.0
n_ldct=39
iwindow=1L
DEVICE, RETAIN=2, DECOMPOSED=0
WINDOW,iwindow,XSIZE=500*fac_window,YSIZE=500*fac_window
loadct,n_ldct

X0=0.10
X1=0.90
Y0=0.10
Y1=0.79
if ( sw_plot_grid ge 1 ) then begin
loadct,0

if ( sw_plot_grid eq 1 ) then begin
plot,plot_xx,plot_yy $
,xrange=[X_min,X_max], xstyle=1 $
,yrange=[Y_min,Y_max], ystyle=1  $
,XTITLE=X_TITLE,YTITLE=Y_TITLE $
,TITLE=VarTitle[VarType]+unit[VarType]+'  ht='+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+'km  UT[hr]='+STRTRIM( string(ut_hr, FORMAT='(F6.2)'),1 )+'.'+TEST $
,POSITION=[X0,Y0,X1,Y1] $
,COLOR=text_color $
,charsize=char_size,charthick=char_thick $
,PSYM=3, SYMSIZE=5.5 ;$
title_plr='ret'

endif  else if ( sw_plot_grid eq 2 ) then begin
;      plot_yy[mp,lp] = mlat_deg[i]
LPmax=30
rmax=90.- plot_yy[0,LPmax-1]
print, 'rmax',rmax
r=fltarr(NMP,LPmax)
theta=fltarr(NMP,LPmax)
for mp=0,NMP-1 do begin
for lp=0,LPmax-1 do begin
r[mp,lp] = 90.- plot_yy[mp,lp] ;degrees from the North mag pole
theta[mp,lp] =  plot_xx[mp,lp] * !PI /180.
if ( r[mp,lp] gt rmax ) then  begin
  r[mp,lp]=0.
  theta[mp,lp]=0.
endif
endfor
endfor
plot,/POLAR, r, theta $
,XSTY=4, YSTY=4 $
;,xrange=[X_min,X_max], xstyle=1 $
;,yrange=[Y_min,Y_max], ystyle=1  $
;,XTITLE=X_TITLE,YTITLE=Y_TITLE $
,TITLE='POLAR PLOT in MAGNETIC COORDINATE' $ ;VarTitle[VarType]+unit+'  ht='+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+'km  UT[hr]='+STRTRIM( string(ut_hr, FORMAT='(F6.2)'),1 ) $
,POSITION=[X0,Y0,X1,Y1] $
,COLOR=text_color $
;,charsize=char_size,charthick=char_thick $
,PSYM=3, SYMSIZE=5.5 ;$
AXIS, 0,0,XAX=0 $
,Xrange=[-rmax,+rmax], Xstyle=1
AXIS, 0,0,YAX=0 $
,Yrange=[-rmax,+rmax], Ystyle=1
title_plr='plr'

endif  ;else if ( sw_plot_grid eq 2 ) then begin

if ( sw_output2file eq 1 ) then begin

   if ( sw_frame eq 0 ) then $  ;magnetic
      title_frame='mag' $
   else if ( sw_frame eq 1 ) then $ ;geographic
      title_frame='geo'

Filename_png=plot_DIR+'ipe_grid_ht'+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+title_plr+'.'+title_frame+'.png'
print,'output to file=',filename_png
output_png, Filename_png
endif

LOADCT, n_ldct
STOP
RETURN
endif ;( sw_plot_grid eq 1 ) then begin

;debug
;mp=25L
;for l=0,ny_max-1 do print,l,plot_zz(mp,l),plot_xx(mp,l),plot_yy(mp,l)
;l=100
;print,'lp==',l
;for mp=0,nx_max-1  do print,mp,plot_zz(mp,l),plot_xx(mp,l),plot_yy(mp,l)

contour,plot_zz,plot_xx,plot_yy $
,/irregular $
,/fill $
,levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
,xrange=[X_min,X_max], /xstyle $
,yrange=[Y_min,Y_max], /ystyle $
,XTITLE=X_TITLE,YTITLE=Y_TITLE $
,TITLE=VarTitle[VarType]+unit[VarType]+'  ht='+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+'km  UT[hr]='+STRTRIM( string(ut_hr, FORMAT='(F6.2)'),1 )+'_'+TEST $
,POSITION=[X0,Y0,X1,Y1] $
,COLOR=text_color $
,charsize=char_size,charthick=char_thick $
,MAX_VALUE= MAX_xymin 

if ( sw_debug eq 1 ) then  print,'MAX=',MAX(plot_zz),' MIN=',MIN(plot_zz)
; add MIN & MAX values
xyouts, 0.6, 0.84 $
, 'MIN='+STRTRIM(STRING( MIN(plot_zz), FORMAT='(E11.3)'),1)+' MAX='+STRTRIM(STRING( MAX(plot_zz), FORMAT='(E11.3)'),1)  $
, charsize=1.0, charthick=1.0, /norm, /noclip

charsize_colorbar=1.5
format_colorbar='(E9.1)'
font=1 ;true-type 
position=[0.10, 0.95, 0.90, 0.98] ;for horizontal bar
COLORBAR, BOTTOM=bottom, CHARSIZE=charsize_colorbar, COLOR=color, DIVISIONS=divisions $
        , FORMAT=format_colorbar, POSITION=position, MAXRANGE=zmax,MINRANGE=zmin $
        , NCOLORS=ncolors,TITLE=title,VERTICAL=vertical,TOP=top,RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
        , _EXTRA=extra, INVERTCOLORS=invertcolors, TICKNAMES=ticknames

if ( sw_output2file eq 1 ) then begin
   if ( sw_frame eq 0 ) then $  ;magnetic
      title_frame='mag' $
   else if ( sw_frame eq 1 ) then $ ;geographic
      title_frame='geo'
Filename_png=plot_DIR+VarTitle[VarType]+'_ht'+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+'_ut'+STRTRIM( string(ut_hr, FORMAT='(F6.2)'),1 )+title_frame+'.png'
output_png, Filename_png
endif


endfor ;VarType=1,1 do begin
end ;pro contour_lon_lat
