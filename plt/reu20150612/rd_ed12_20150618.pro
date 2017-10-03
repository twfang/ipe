;date: 20150618
;version: 4: plot input
pro rd_ed12_20150618

runid='ipe_S_13871'

; path to the output files
input_DIR0='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r311.1tmp/trunk/run/'+runid

nmp=80L
nlp=170L
mlon90_2d=fltarr(nmp,nlp*2)
mlat90_2d=fltarr(nmp,nlp*2)
;calculate mlat90_2d
mlat90_1=fltarr(nlp*2)
title_res='low'
sw_debug=1L
get_gl, mlat90_1,title_res,sw_debug
print, 'mlat90_1', mlat90_1
for mp=0,nmp-1 do begin
mlat90_2d[mp,0:nlp*2-1]=mlat90_1[0:nlp*2-1]
endfor

;calculate mlon90_2d
dlonm90km=4.50 ;deg
mlon90=findgen(nmp)*dlonm90km
print,'mlon90',mlon90
for lp=0,nlp*2-1 do begin
   for mp=0,nmp-1 do begin
      mlon90_2d[mp,lp]=mlon90[mp]
   endfor ;mp
endfor ;lp

input_DIR='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/reu/tmp20130703reu/trunk/run/ipe_640_18702'

openr, LUN6, input_DIR+'/ut_rec', /GET_LUN ;utime
openr, LUN8, input_DIR+'/fort.2008', /GET_LUN ;ed1_90
openr, LUN9, input_DIR+'/fort.2009', /GET_LUN ;ed2_90


formatE='(20E12.4)'
ed190=fltarr(nmp,nlp*2) ;[mV/m]
ed290=fltarr(nmp,nlp*2) ;[mV/m]


nmin=1-1L
nmax=97-1L;97L
vExBUp=fltarr(nmax+1)
vExBe =fltarr(nmax+1)
lthr=fltarr(nmax+1)
for n_read=nmin,nmax do begin

;read ed190
   edum=fltarr(2,nlp,nmp)
   readf, LUN8, edum,  FORMAT=formatE ;ed190
   for mp=0,nmp-1 do begin
      ii=-1
      for lp=0,nlp-1 do begin
         ii=ii+1
         ed190[mp,ii] = edum[1,lp,mp] ;SH
      endfor                          ;lp
      for lp=nlp-1,0,-1 do begin
         ii=ii+1
         ed190[mp,ii] = edum[0,lp,mp] ;NH
      endfor                          ;lp
   endfor                             ;mp

;read ed290
   edum=fltarr(2,nlp,nmp)
   readf, LUN9, edum,  FORMAT=formatE ;ed290
   for mp=0,nmp-1 do begin
      ii=-1
      for lp=0,nlp-1 do begin
         ii=ii+1
         ed290[mp,ii] = edum[1,lp,mp]*(-1.) ;SH
      endfor                          ;lp
      for lp=nlp-1,0,-1 do begin
         ii=ii+1
         ed290[mp,ii] = edum[0,lp,mp]*(-1.) ;NH
      endfor                          ;lp
   endfor                             ;mp
 


   v_e1 =   ed290[mp0,lp0]*1.E-3 / Be3 ;!(4.18) +mag-east 
   v_e2 = - ed190[mp0,lp0]*1.E-3 / Be3 ;!(4.19) +down/equatorward
   
   vExBUp[n_read] = v_e2 * (-1.)  ;down-->upward
   vExBe[ n_read] = v_e1          ;
   print, 'n_read=', n_read,' ed1[mV/m]=',ed190[mp0,lp0],' ed2[mV/m]=',ed290[mp0,lp0]
   print,' vEXBup[m/s]=', vExBup[n_read],' vEXBe[m/s]=', vExBe[n_read]


;read Universal Time (UT)
rec=0L
utime=0L
ut0=168.
   readf, LUN6, rec, utime           ;in seconds
   uthr = utime / 3600. - ut0
;calculate Local Time (LT)
   lthr[n_read] = uthr - 5.
   print, 'utime=',utime,' UT hr=', uthr,' LT=',lthr[n_read]
   
endfor                          ;n_read=nmin,nmax

;close file
free_lun,lun8


;plot EXB drift as a function of LT
;set up a display window
DEVICE, RETAIN=2, DECOMPOSED=0
WINDOW,1,XSIZE=1000,YSIZE=1000
loadct, 0  ;black and white
!p.multi=[0,3,2,0]

;plot vEXB upward component
plot,lthr,vExBup $
, yrange=[-90. , +90. ],  ystyle=1  $
, xrange=[ MIN(lthr), MIN(lthr) ],  xstyle=1  $
,XTITLE = 'LT[hr] ', YTITLE = 'VEXB[m/s]' $
, charsize =2.5, charthick=2.5 $
, TITLE=' ' $
, linestyle=0 $
, thick=2.0 $
, color=250.

;overplot the eastward drift 
oplot,lthr,vExBe $
, linestyle=2 $
, thick=2.0 $
, color=250.



n_levels=100
zmax=MAX(ed190)
zmin=MIN(ed190)
X_min=MIN(mlon90_2d)
X_max=MAX(mlon90_2d)
Y_min=MIN(mlat90_2d)
Y_max=MAX(mlat90_2d)
text_color=160
char_size =2.0
char_thick=1.5
loadct,39
; plot ed1 in color contour
contour,z_data,mlon90_2d,mlat90_2d $
,/irregular $
,/fill $
, levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
, xrange=[X_min,X_max], /xstyle  $
, yrange=[Y_min,Y_max], /ystyle  $
,XTITLE = 'mlon90', YTITLE = 'mlat90' $ 
,TITLE = 'Ed1 ' $
, COLOR=text_color $
, charsize = char_size, charthick = char_thick


; Add a colorbar
charsize_colorbar=3.0
format_colorbar='(E10.2)' ;No+
font=1 ;true-type 
position=[0.38, 0.498, 0.95, 0.503]  ;for horizontal bar
COLORBAR, BOTTOM=bottom, CHARSIZE=charsize_colorbar, COLOR=color, DIVISIONS=divisions $
        , FORMAT=format_colorbar, POSITION=position, MAXRANGE=zmax,MINRANGE=zmin $
        , NCOLORS=ncolors,TITLE=title,VERTICAL=vertical,TOP=top,RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
        , _EXTRA=extra, INVERTCOLORS=invertcolors, TICKNAMES=ticknames



;save figure in a file
output_png, 'vexb_ed1.png'

openw,lun5,'output_lthr', /GET_LUN
printf, lun5,lthr
free_lun,lun5

openw,lun5,'output_exbup', /GET_LUN
printf, lun5,vExBup
free_lun,lun5

openw,lun5,'output_exbe', /GET_LUN
printf, lun5,vExBe
free_lun,lun5

print, 'program rd_ed12 finished!'
end ;pro rd_ed12
