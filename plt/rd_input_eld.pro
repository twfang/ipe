;date: 20150618
;version: 4: plot input to the electrodynamo solver
pro rd_input_20150618

runid='ipe_S_13871'

; path to the output files
input_DIR='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r311.1tmp/trunk/run/'+runid

kmlonp1=81L
kmlat  =97L
formatE='(20E12.4)'

;read mlat
mlat=fltarr(kmlat)
openr, LUN8, input_DIR+'/fort.4020', /GET_LUN ;mlat
readf, LUN8, mlat,  FORMAT=formatE
print,'mlat',mlat
;read mlon
mlon=fltarr(kmlonp1)
readf, LUN8, mlon,  FORMAT=formatE
print,'mlon',mlon
;close file
free_lun,lun8

;open time dependent output files
openr, LUN7, input_DIR+'/ut_rec', /GET_LUN ;utime
openr, LUN8, input_DIR+'/fort.4022', /GET_LUN ;zigm11
openr, LUN9, input_DIR+'/fort.4023', /GET_LUN ;zigm22
openr, LUN10, input_DIR+'/fort.4024', /GET_LUN ;zigmc
openr, LUN11, input_DIR+'/fort.4025', /GET_LUN ;zigm2
openr, LUN12, input_DIR+'/fort.4026', /GET_LUN ;rim


nmin=1-1L
nmax=1-1L;97L
for n_read=nmin,nmax do begin
   print, 'n_read=', n_read

;read zigm11
   zigm11=fltarr(kmlonp1,kmlat)
   readf, LUN8, zigm11,  FORMAT=formatE ;zigm11
;...
;...

; read rim
   rim=fltarr(kmlonp1,kmlat,2)

;read Universal Time (UT)
rec=0L
utime=0L
   readf, LUN7, rec, utime           ;in seconds
   uthr = utime / 3600. 
   print, 'utime=',utime,' UT hr=', uthr  

;choose zigm11 for plotting
z_data=zigm11
var_title='zigm11'

;set up a display window
DEVICE, RETAIN=2, DECOMPOSED=0
WINDOW,1,XSIZE=1000,YSIZE=800
loadct, 39  ;black and white
!p.multi=[0,3,2,0]

;set up parameters required for contour plot
n_levels=100
zmax=MAX(z_data)
zmin=MIN(z_data)
X_min=MIN(mlon)
X_max=MAX(mlon)
Y_min=MIN(mlat)
Y_max=MAX(mlat)
text_color=255.;160
char_size =2.0
char_thick=1.5
loadct,39

; plot input(fli) in color contour
contour,z_data,mlon,mlat $
,/fill $
, levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
, xrange=[X_min,X_max], /xstyle  $
, yrange=[Y_min,Y_max], /ystyle  $
,XTITLE = 'mlon[deg]', YTITLE = 'mlat[deg]' $ 
,TITLE = var_title $
, COLOR=text_color $
, charsize = char_size, charthick = char_thick 

; Add a colorbar
charsize_colorbar=3.
format_colorbar='(E8.0)' ;No+
font=1 ;true-type 
x0=0.025
dx=0.295
y0=0.5
dy=0.005
position=[x0, y0, (x0+dx), (y0+dy)]  ;for horizontal bar
COLORBAR, BOTTOM=bottom, CHARSIZE=charsize_colorbar, COLOR=color, DIVISIONS=divisions $
        , FORMAT=format_colorbar, POSITION=position, MAXRANGE=zmax,MINRANGE=zmin $
        , NCOLORS=ncolors,TITLE=title,VERTICAL=vertical,TOP=top,RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
        , _EXTRA=extra, INVERTCOLORS=invertcolors, TICKNAMES=ticknames


; print out UT 
xyouts, 0.010, 0.980 $
, 'UT[hrs]='+STRTRIM( string(uthr, FORMAT='(F5.1)'), 1)  $
, charsize =1.5, charthick=1.0 $
, /norm, /noclip

; print out runid
xyouts, 0.8, 0.07 $
, runid  $
, charsize =1.5, charthick=1.0 $
, /norm, /noclip


;save figure in a file
output_png, 'input_ut'+STRTRIM( string(uthr, FORMAT='(F5.1)'), 1)+'.png'

endfor                          ;n_read=nmin,nmax

;close file
free_lun,lun8
free_lun,lun9
free_lun,lun10
free_lun,lun11
free_lun,lun12



print, 'program rd_input finished!'
end ;pro rd_input
