pro plt_r_sza_ne $
  , JMIN_IN,JMAX_IS,Z_km, sza_rad $
  , xionn_m3, te_ti_k $
  , sw_debug,mlat_deg,ut_hr,input_DIR0,plot_DIR

swParam=1L ;0=Ne;1=Te
sw_save=1L
sw_plot_contour=1L
utHrPlt0=getenv('utHrPlt0')
jthMax=getenv('jthMax')
sw_output2file=getenv('sw_output2file') 

;dbg output file
;before SH
if sw_debug eq 1 then begin
   lunTmp0=2000L
   openw,luntmp0,'tmp0.dat', /GET_LUN
;after SH just before save
   lunTmp1=2001L
   openw,luntmp1,'tmp1.dat', /GET_LUN
endif ;sw_debug eq 1 then begin







if sw_debug eq 1 then  print,'!dbg plt_r_sza', sw_hem
size_result = SIZE(JMIN_IN)
NLP=size_result[1]
if ( sw_debug eq 1 ) then  print, 'plt_r_sza:NLP=', NLP

size_result = SIZE(xionn_m3)
NMP=size_result[3]
if ( sw_debug eq 1 ) then  print, 'plt_r_sza:NMP=', NMP
ISPEC=size_result[1]
if ( sw_debug eq 1 ) then  print, 'plt_r_sza:ISPEC=', ISPEC


Re_km=6.3712E+03;km
fac_den=getenv('factor_density') 

d_sza=0.1
sza_min=70.
sza_max=130.
imax=FIX((sza_max-sza_min)/d_sza)+1
r_min=1.
r_max=2.256 
jmax=202L
y=fltarr(jmax)
if sw_debug eq 1 then  print, 'imax', imax,'jmax',jmax
missingValue=-9999.99
nhem=2L
 z_sav=fltarr(nhem,imax,jmax)
nz_sav=LONarr(nhem,imax,jmax)
lp_max=22L ;mlat=64.7976deg, assuming inside polar cap


for sw_hem=0,nhem-1 do begin ;0='NH';1='SH'
   for mp=0,nmp-1 do begin

   sza_degMax=-9999.99
   mlat_degMax=-9999.99
   z_kmMax=-9999.99
   for lp=0,lp_max do begin

;i: sza at the bottom of the flux tube
      if sw_hem eq 0 then $     ;MH
         ipts=jmin_in[lp]-1L $
      else if sw_hem eq 1 then $ ;SH
         ipts=jmax_is[lp]-1L

      sza_deg = sza_rad[ipts,mp]*180./!pi ;sza [rad]-->deg
      if sza_deg lt sza_min or sza_deg gt sza_max then continue
      i=FIX((sza_deg-sza_min)/d_sza)
      if sw_debug eq 1 then print,'i=', i, sza_deg,lp,mlat_deg[ipts],mp,z_km[ipts]
      if sza_deg gt sza_degMax then begin
         sza_degMax=sza_deg
         mlat_degMax=mlat_deg[ipts]
         z_kmMax=z_km[ipts]
      endif



midpoint =jmin_in[lp]+ (jmax_is[lp]-jmin_in[lp]    )/2  -1 
if sw_hem eq 0L then begin
;NH
   iptsMin=jmin_in[lp]-1L
   iptsMax=jmin_in[lp]-1L+jmax
   istep=+1L
endif else if sw_hem eq 1L then begin
;SH
   iptsMin=jmax_is[lp]-1L
   iptsMax=jmax_is[lp]-jmax
   istep=-1L
endif

      for ipts=iptsMin,iptsMax,istep do begin

         if sw_hem eq 0L then $
            j=ipts-iptsMin $
         else if sw_hem eq 1L then $
            j=-ipts+iptsMin

        r=(z_km[ipts]+re_km)/re_km ;r
        if r gt r_max then continue
        y[j]=r
        if sw_debug eq 1 then   print,sw_hem,'      j=', j,r,ipts,z_km[ipts]

        if swParam eq 0L then begin ;0=Ne;1=Te
;derive Ne
           Nel=0.00000
           for jth=0,jthMax do begin
              Nel = Nel + xionn_m3[jth,ipts,mp]
           endfor

           if ( Nel gt 0. ) then $
              alog10Nel = ALOG10( Nel*fac_den ) $ ;ne [log10 cm-3]
           else begin
              print,'!INVALID Ne! SKIP',nel,ipts,lp,mp
                                ;alog10Nel = -6.
              CONTINUE          ;goto next ipts
           endelse
        endif else if swParam eq 1L then begin ;1=Te
           alog10Nel=te_ti_k[3-1,ipts,mp]
        endif                   ;swParam eq 0L then begin ;0=Ne;1=Te

                                ;dbg z_sav [i,j]= alog10Nel
        tmp=z_sav[sw_hem,i,j]
        z_sav [sw_hem,i,j] = tmp + alog10Nel
        tmp=nz_sav[sw_hem,i,j]  
        nz_sav[sw_hem,i,j] = tmp + 1
        if sw_debug eq 1 then print,i,j,z_sav[sw_hem,i,j],nz_sav[sw_hem,i,j],ipts,lp,mp

     endfor ;ipts

   endfor                          ;lp
   if sw_debug eq 1 then print,sw_hem,mp,'sza Max=',sza_degMax,mlat_degMax,z_kmMax
endfor                             ;mp

;take the average
for j=0,jmax-1 do begin
   for i=0,imax-1 do begin
      tmp=z_sav[sw_hem,i,j]
      if nz_sav[sw_hem,i,j] gt 0 then $
        z_sav[sw_hem,i,j] = tmp / FLOAT(nz_sav[sw_hem,i,j]) $
      else $
        z_sav[sw_hem,i,j] = missingValue
      if sw_debug eq 1 then print,i,j,z_sav[sw_hem,i,j],tmp,nz_sav[sw_hem,i,j]
   endfor                     ;i
endfor                          ;j


if sw_plot_contour eq 1 then begin

z=fltarr(imax,jmax)
for j=0,jmax-1 do begin
   for i=0,imax-1 do begin
      z[i,j]=z_sav[sw_hem,i,j]
   endfor                       ;i
endfor                          ;j
print,'sw_hem',sw_hem,'  MAX Z=',MAX(z),MIN(z)
text_color=255.
char_size=3.0
char_thick=3.0
n_ldct=39

	iwindow=1L
	DEVICE, RETAIN=2, DECOMPOSED=0
	WINDOW,iwindow,XSIZE=1000,YSIZE=1000
	!p.multi=[0,1,1,0]
	loadct,n_ldct

x=findgen(imax)*d_sza + sza_min
if sw_debug eq 1 then  print,' x=', x
if sw_debug eq 1 then  print,' y=', y

if swParam eq 0 then begin
  zmax=5.5
  zmin=0.0
endif else if swParam eq 1 then begin
  zmax=12600.
  zmin=180.0
endif
 
n_levels=100
x_min=sza_min
x_max=sza_max
y_min=r_min
y_max=r_max
x_title='SZA [deg]'
y_title='R'
pltUt=ut_hr-utHrPlt0
if pltUt ge 10. then $
  titleUt=STRTRIM( string(pltUt, FORMAT='(F6.2)'),1 ) $
else $  ;if pltUt ge 10. then $
  titleUt='0'+STRTRIM( string(pltUt, FORMAT='(F6.2)'),1 )


if sw_hem eq 0L then $
   titleHem='NH' $
else if sw_hem eq 1L then $
   titleHem='SH' 

Y0=0.12
Y1=0.9

;dbg20161206
if sw_debug eq 1 and sw_hem eq 0 then begin
   for i=0,imax-1 do begin 
      for j=0,jmax-1 do begin 
;      printf,lunTmp0,i,j,z_sav[0,i,j]
         printf,lunTmp0,i,j,z[i,j]
      endfor                    ;j
   endfor                       ;i
   free_lun,luntmp0
   for i=0,1 do print,'0ihem=',i,' Z MAX',max(z_sav[i,*,*],Max_Subscript),min(z_sav[i,*,*]), Max_Subscript
   print,sw_hem,' x',size(x),MAX(x), x
   print,sw_hem,' y',size(y),MAX(y), y
endif                           ;sw_debug eq 1


if swParam eq 0 then $
   titleVar = 'Ne [log!D10!N cm!U-3!N]' $
else if swParam eq 1 then $
   titleVar = 'Te [K]'

;contour plot fig3
contour,z,x,y $
,/fill $
,levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
,xrange=[X_min,X_max], /xstyle $
,yrange=[Y_min,Y_max], /ystyle $
,XTITLE=X_TITLE,YTITLE=Y_TITLE $
,TITLE=titleVar+': '+titleHem+' UT[hr]='+titleUt $
;TEST+'_'+rundir $
,POSITION=[0.125,Y0,0.9,Y1] $;X0,Y0,X1,Y1] $
,COLOR=text_color $
,charsize=char_size,charthick=char_thick ;$
;,MIN_VALUE= (missingValue+0.1)



xyouts, 0.38, 0.010 $
,input_DIR0 $
,charsize=0.895, charthick=0.8, /norm, /noclip


title = ' ';colorbar title'
font = 1  ;True-Type: 1.
if swParam eq 0 then $
   format = '(f7.2)' $
else if swParam eq 1 then $
   format = '(f7.0)'


vertical = 'vertical'
position = [0.97, Y0, 0.98, Y1] ;for a vertical bar 

; add colorbar
COLORBAR, BOTTOM=bottom, CHARSIZE=char_size, COLOR=color, DIVISIONS=divisions $
        , FORMAT=format, POSITION=position, MAXRANGE=zmax, MINRANGE=zmin $
        , NCOLORS=ncolors, TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
        , _EXTRA=extra, INVERTCOLORS=invertcolors,  TICKNAMES=ticknames

filename_image=plot_DIR+'sza'+titleHem+'_ut'+titleUt+'p'+STRTRIM( string(swParam, FORMAT='(i1)'),1 )+'.png' 
if sw_output2file eq 1 then  output_png, filename_image 
endif ;sw_plot_contour eq 1 then 
endfor ;sw_hem


;dbg20161206
if sw_debug eq 1 then begin
   for i=0,imax-1 do begin 
      for j=0,jmax-1 do begin 
         printf,lunTmp1,i,j,z_sav[0,i,j]
      endfor                    ;j
   endfor                       ;i
   free_lun,luntmp1
;print,'imax',imax,'jmax',jmax
   for i=0,1 do print,'1ihem=',i,' Z MAX',max(z_sav[i,*,*],Max_Subscript),min(z_sav[i,*,*]), Max_Subscript
   print,sw_hem,' 1x',size(x),MAX(x), x
   print,sw_hem,' 1y',size(y),MAX(y), y
   if sw_hem eq 0 then begin
      yNH=y
   endif else if sw_hem eq 1 then begin
      for j=0,jmax-1 do begin
         dif=(y[j]-yNH[j])
         print,'check dif: j=', j, dif
         if dif ne 0.0 then STOP
      endfor                    ;j
   endif ;sw_hem
endif    ;sw_debug eq 1 then 


;save
if sw_save eq 1 then begin
   flnm_sav=plot_DIR+'sza'+titleUt+'p'+STRTRIM( string(swParam, FORMAT='(i1)'),1 )+'.sav' 
   print,'saving file=', flnm_sav
   save, z_sav,x,y,missingValue, /variables, filename=flnm_sav
endif ;sw_save


if sw_debug eq 1 then print, "plt_r_sza_ne successfully finished!"
end ;pro plt_r_sza_ne
