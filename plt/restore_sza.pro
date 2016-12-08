pro restore_sza
sw_plot_contour=1L;1
sw_debug=0L
plot_DIR='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/fig/20161128kitamura/1480638394_ipe_theia_intel_parallel2_80/'
n_read_max=1L;97L
dUt=0.25000
ut0=0.00

nhem=2L
imax=601L
jmax=202L
print,'0imax',imax,'jmax',jmax
zAl=fltarr( nhem,imax,jmax)
nzAl=LONarr(nhem,imax,jmax)

for nt=0,n_read_max-1 do begin

   pltUt=ut0+dUt*FLOAT(nt)

   if pltUt ge 10. then $
      titleUt=STRTRIM( string(pltUt, FORMAT='(F6.2)'),1 ) $
   else $                       ;if pltUt ge 10. then $
      titleUt='0'+STRTRIM( string(pltUt, FORMAT='(F6.2)'),1 )

;restore
   flnm_sav=plot_DIR+'sza'+titleUt+'.sav' 
   print,'nt=',nt,pltUt,' restore file=', flnm_sav
   restore, flnm_sav
;dbg NH
print,'1imax',imax,'jmax',jmax
for i=0,1 do print,'ihem=',i,' Z MAX',max(z_sav[i,*,*],Max_Subscript),min(z_sav[i,*,*]), Max_Subscript
lunTmp2=2002L
openw,luntmp2,'tmp2.dat', /GET_LUN
for i=0,imax-1 do begin 
   for j=0,jmax-1 do begin 
      printf,lunTmp2,i,j,z_sav[0,i,j]
   endfor                       ;j
endfor  
free_lun,luntmp2

size_result=SIZE(z_sav)
nhem=size_result[1]
imax=size_result[2]
jmax=size_result[3]
if sw_debug eq 1 then print,nhem,imax,jmax


for ihem=0,nhem-1 do begin
   for i=0,imax-1 do begin
      for j=0,jmax-1 do begin



         tmp =zAl[ ihem,i,j]
         ntmp=nzAl[ihem,i,j]

         if z_sav[ihem,i,j] gt missingValue then begin
            zAl[ ihem,i,j]= tmp+z_sav[ihem,i,j] 
            nzAl[ihem,i,j]=ntmp+1            
if sw_debug eq 1 AND ihem eq 0 then  print,'0ihem=',ihem,i,j, zAl[ ihem,i,j],tmp,z_sav[ihem,i,j],nzAl[ihem,i,j]
         endif 

      endfor                    ;j
   endfor                       ;i
endfor                          ;ihem

;time loop end
endfor                          ;nt
;STOP ;debug

;average
for ihem=0,nhem-1 do begin
   for i=0,imax-1 do begin
      for j=0,jmax-1 do begin
         tmp=zAl[ihem,i,j]
         if nzAl[ihem,i,j] gt 0 then $
            zAl[ihem,i,j] = tmp / FLOAT(nzAl[ihem,i,j]) $
         else $
            zAl[ihem,i,j] = missingValue
;         if sw_debug eq 1 then $
         if sw_debug eq 1 AND ihem eq 0 then $
print,'1ihem=',ihem,i,j,zAl[ihem,i,j],tmp,nzAl[ihem,i,j]
      endfor                    ;j
   endfor                       ;i
endfor                          ;ihem

if sw_plot_contour eq 1 then begin
;plot total Z
   for ihem=0,nhem-1 do begin

   z=fltarr(imax,jmax)
   for i=0,imax-1 do begin
      for j=0,jmax-1 do begin
;         z[i,j]=zAl[ihem,i,j]
         z[i,j]=z_sav[ihem,i,j] ;;dbg!!!
      endfor                    ;j
   endfor                       ;i
   print,ihem,'MAX Z=',MAX(z),MIN(z)

;contour plot
text_color=255.
char_size=3.0
char_thick=3.0
n_ldct=39

	iwindow=(ihem+2L)
	DEVICE, RETAIN=2, DECOMPOSED=0
	WINDOW,iwindow,XSIZE=1000,YSIZE=1000
	!p.multi=[0,1,1,0]
	loadct,n_ldct

;x=findgen(imax)*d_sza + sza_min
if sw_debug eq 1 then  print,' x=', x
if sw_debug eq 1 then  print,' y=', y

zmax=5.5
zmin=0.0
n_levels=100
x_min=70.;sza_min
x_max=130.;sza_max
y_min=1.;r_min
y_max=2.256;r_max
x_title='SZA [deg]'
y_title='R'


if ihem eq 0 then $
   titleHem='NH' $
else if ihem eq 1 then $
   titleHem='SH' 

Y0=0.12
Y1=0.9

;dbg NH
if ihem eq 0 then begin 
lunTmp3=2003L
openw,luntmp3,'tmp3.dat', /GET_LUN
for i=0,imax-1 do begin 
   for j=0,jmax-1 do begin 
      printf,lunTmp3,i,j,z[i,j]
   endfor                       ;j
endfor  
free_lun,luntmp3
endif ;ihem
for i=0,1 do print,'ihem=',i,' Z MAX',max(z_sav[i,*,*],Max_Subscript),min(z_sav[i,*,*]), Max_Subscript


print,ihem,' 1x',size(x),MAX(x), x
print,ihem,' 1y',size(y),MAX(y), y
;contour plot fig3
contour,z,x,y $
,/fill $
,levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
,xrange=[X_min,X_max], /xstyle $
,yrange=[Y_min,Y_max], /ystyle $
,XTITLE=X_TITLE,YTITLE=Y_TITLE $
,TITLE='Ne [log!D10!N cm!U-3!N]: '+titleHem $
,POSITION=[0.125,Y0,0.9,Y1] $;X0,Y0,X1,Y1] $
,COLOR=text_color $
,charsize=char_size,charthick=char_thick ;$
;,MIN_VALUE= (missingValue+0.1)



;xyouts, 0.38, 0.010 $
;,input_DIR0 $
;,charsize=0.895, charthick=0.8, /norm, /noclip


title = ' ';colorbar title'
font = 1  ;True-Type: 1.
format = '(f7.2)'

vertical = 'vertical'
position = [0.97, Y0, 0.98, Y1] ;for a vertical bar 

; add colorbar
COLORBAR, BOTTOM=bottom, CHARSIZE=char_size, COLOR=color, DIVISIONS=divisions $
        , FORMAT=format, POSITION=position, MAXRANGE=zmax, MINRANGE=zmin $
        , NCOLORS=ncolors, TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
        , _EXTRA=extra, INVERTCOLORS=invertcolors,  TICKNAMES=ticknames

filename_image=plot_DIR+'sza'+titleHem+'Total.png' 
output_png, filename_image


endfor;ihem

endif ;sw_plot_contour eq 1 then begin
 print,' end pro restore_sza'
  end ;pro restore_sza
