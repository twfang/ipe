;plotting does not work!!!
;20161213 special version to plot Te and Ti at the same time
;by restoring the saved Total Ti data
pro restore_sza_tcombined


;restore grid for z_km
filename_grid_sav='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/grid/plt/plasma_grid.2xdyn.sav'
print,'restoring the grid file=',filename_grid_sav
restore, filename=filename_grid_sav
if sw_debug eq 1 then print, jmin_in[lp0-1],jmax_is[lp0-1]
for i=0,202 do print,i, z_km[i]
;--

fac_window=1.0
swParam=2L   ;0:Ne;1:Te;2:Ti
sw_saveTotal=0L
sw_plot_contour=1L;1
sw_debug=0L
rundir=$
'1481593483_ipe_theia_intel_parallel2_80' ;v31new run with no conjugate PE
;'1481414104_ipe_theia_intel_parallel2_80' ;v31new run with sw_perp_trans=0 output/2hr
;'1481413483_ipe_theia_intel_parallel2_80' ;v30new run with sw_perp_trans=0
;'1480638394_ipe_theia_intel_parallel2_80' ;v28old run with sw_perp_trans=1
plot_DIR0= $
'/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/fig/20161128kitamura/'+rundir+'/'
;;theia
;'/home/Naomi.Maruyama/p2/ipe/20161128kitamura/Ne/'+rundir+'/' ;jet
plot_DIR1=plot_DIR0+'total/'
n_read_max=$
;1L ;debug
13L;new v31
;25L;new v30
;97L;old v28
dUt=$
2.0;new v31
;1.0;new v30
;0.25000;old v28
ut0=0.00

nhem=2L
imax=601L
jmax=202L
if sw_debug eq 1 then  print,'0imax=',imax,' jmax=',jmax
zAl=fltarr( nhem,imax,jmax)
nzAl=LONarr(nhem,imax,jmax)



for nt=0,n_read_max-1 do begin

   pltUt=ut0+dUt*FLOAT(nt)

   if pltUt ge 10. then $
      titleUt=STRTRIM( string(pltUt, FORMAT='(F6.2)'),1 ) $
   else $                       ;if pltUt ge 10. then $
      titleUt='0'+STRTRIM( string(pltUt, FORMAT='(F6.2)'),1 )

   ;if swParam eq 0 then $
   ;   titlePar='' $
   ;else if swParam ge 1 then $
      titlePar='p'+STRTRIM( string(swParam, FORMAT='(i1)'),1 )

;restore
   flnm_sav=plot_DIR0+'sza'+titleUt+titlePar+'.sav' 
   print,'nt=',nt,pltUt,' restore file=', flnm_sav
   restore, flnm_sav
;dbg NH
if sw_debug eq 1 then begin
   print,'1imax',imax,' jmax',jmax
   for i=0,1 do print,'ihem=',i,' Z MAX',max(z_sav[i,*,*],Max_Subscript),min(z_sav[i,*,*]), Max_Subscript
   lunTmp2=2002L
   openw,luntmp2,'tmp2.dat', /GET_LUN
   for i=0,imax-1 do begin 
      for j=0,jmax-1 do begin 
         printf,lunTmp2,i,j,z_sav[0,i,j]
      endfor                    ;j
   endfor  
   free_lun,luntmp2
endif ;sw_debug

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
;note: ihem=2: average both hems
;   for ihem=0,nhem do begin
   for ihem=nhem,nhem do begin ;plot only NH+SH combined data

   z=fltarr(imax,jmax)
   for i=0,imax-1 do begin
      for j=0,jmax-1 do begin

         if ihem lt nhem then $
           z[i,j]=zAl[ihem,i,j] $
         else if ihem eq nhem then $
           z[i,j]=TOTAL(zAl[*,i,j])*0.50

      endfor                     ;j
   endfor                       ;i
   print,'ihem=',ihem,' MAX Z=',MAX(z),MIN(z)

;contour plot
text_color=255.
char_size=3.0
char_thick=3.0
n_ldct=39

	iwindow=ihem
	DEVICE, RETAIN=2, DECOMPOSED=0
	WINDOW,iwindow,XSIZE=2000*fac_window,YSIZE=1000*fac_window
	!p.multi=[0,3,1,0]
	loadct,n_ldct

;x=findgen(imax)*d_sza + sza_min
if sw_debug eq 1 then  print,' x=', x
if sw_debug eq 1 then  print,' y=', y


if swParam eq 0 then begin
   zmax=5.5
   zmin=0.0
endif else if swParam ge 1 then begin
   zmax=6000.
   zmin=200.
endif

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
   titleHem='SH' $
else if ihem eq 2 then $
   titleHem='NH+SH' 

Y0=0.12
Y1=0.9

;dbg NH
if sw_debug eq 1 then begin 
   lunTmp3=2003L
   openw,luntmp3,'tmp3.dat', /GET_LUN
   for i=0,imax-1 do begin 
      for j=0,jmax-1 do begin 
         printf,lunTmp3,i,j,z[i,j]
      endfor                    ;j
   endfor  
   free_lun,luntmp3
   for i=0,1 do print,'ihem=',i,' Z MAX',max(z_sav[i,*,*],Max_Subscript),min(z_sav[i,*,*]), Max_Subscript
   print,ihem,' 1x',size(x),MAX(x), x
   print,ihem,' 1y',size(y),MAX(y), y
endif                           ;sw_debug


if swParam eq 0 then begin
  titlePar='Ne' 
  titleUni='log!D10!N cm!U-3!N' 
endif else if swParam eq 1 then begin
  titlePar='Te'
  titleUni='K'
endif else if swParam eq 2 then begin
  titlePar='Ti'
  titleUni='K'
endif

;contour plot fig3
contour,z,x,y $
,/fill $
,levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
,xrange=[X_min,X_max], /xstyle $
,yrange=[Y_min,Y_max], /ystyle $
,XTITLE=X_TITLE,YTITLE=Y_TITLE $
,TITLE=titlePar+' ['+titleUni+']: '+titleHem $
;,POSITION=[0.125,Y0,0.9,Y1] $;X0,Y0,X1,Y1] $
,POSITION=[0.125,Y0,0.49,Y1] $;X0,Y0,X1,Y1] $
,COLOR=text_color $
,charsize=char_size,charthick=char_thick ;$
;,MIN_VALUE= (missingValue+0.1)





title = ' ';colorbar title'
font = 1  ;True-Type: 1.
if swParam eq 0 then $
   format = '(f7.2)' $
else if swParam ge 1 then $
   format = '(f7.0)'


vertical = 'vertical'
position = [0.52, Y0, 0.53, Y1] ;for a vertical bar 

; add colorbar
COLORBAR, BOTTOM=bottom, CHARSIZE=char_size, COLOR=color, DIVISIONS=divisions $
        , FORMAT=format, POSITION=position, MAXRANGE=zmax, MINRANGE=zmin $
        , NCOLORS=ncolors, TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
        , _EXTRA=extra, INVERTCOLORS=invertcolors,  TICKNAMES=ticknames

;add ht profiles
;chose sza
;for i=0,150 do print, i, x[i]
i1=50L
i2=550L
yz=y
for j=0,jmax-1 do yz[j]=z_km[j]
X0=0.58
dX=0.2
X1=X0+dX
x_title=titlePar+' ['+titleUni+']'
y_title='Alt [km]'
print,'i1',x[i1], MAX(z[i1,*]), MIN(z[i1,*])
if swParam eq 0 then begin
  y_max=8000. 
  x_max=6.
  endif else begin ;$                   ;if swParam eq 0 then $
  y_max=3500.
  x_max=6000.
endelse
y_min=200.
x_min=0.
;1)sza=75
plot, z[i1,*],yz $
,xrange=[X_min,X_max], /xstyle, yrange=[Y_min,Y_max], /ystyle $
;,TITLE='SZA[deg]=',+STRTRIM( string(x[i1], FORMAT='(f6.0)'),1 ) $
,XTITLE=X_TITLE $
,YTITLE=Y_TITLE $
,charsize=char_size,charthick=char_thick $
,POSITION=[X0,Y0,X1,Y1]


;20161213 changes from here, ,,,in order to plot Te and ADD Ti at the
;same time 
X01=X0+dX
X11=X0+dX
xi2=x[i2]
z2=fltarr(jmax)
for j=0,jmax-1 do z2[j]=z[i2,j]
print,'i2',xi2, MAX(z2), MIN(z2)

; safe to restore Ti here
restore,'/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/fig/20161128kitamura/1481593483_ipe_theia_intel_parallel2_80/total/sza24.00TiTotal.sav'

;over-plot Ti sza=75
print,i1, z[i1,*],yz
oplot, z[i1,*],yz $
 ,linestyle=2


;2)sza=125
plot, z2,yz $
,xrange=[X_min,X_max], /xstyle, yrange=[Y_min,Y_max], /ystyle $
,charsize=char_size,charthick=char_thick $
,POSITION=[X01,Y0,X11,Y1]


;over-plot Ti sza=125
oplot, z[i2,*],yz $
 ,linestyle=2


xyouts, 0.90, 0.030 $
;, 'MAX:'+STRTRIM( string(MAX(z),  FORMAT=format) ,1)+' MIN:'+STRTRIM( string(MIN(z),  FORMAT=format) ,1)+' ['+titleUni+']' $
, 'MAX: '+STRTRIM( string(MAX(z2),  FORMAT=format) ,1)+' ['+titleUni+']' $
,charsize=1., charthick=1., /norm, /noclip


xyouts, 0.38, 0.010 $
,input_DIR0 $
,charsize=0.895, charthick=0.8, /norm, /noclip


filename_image=plot_DIR1+'sza'+titlePar+titleHem+'TotalTei.png' 
output_png, filename_image


endfor;ihem

endif ;sw_plot_contour eq 1 then begin



;save
if sw_saveTotal eq 1L then begin
   flnm_sav=plot_DIR1+'sza'+titleUt+titlePar+'Total.sav' 
   print,titlePar,' saving Total file=', flnm_sav
   save, zAl,x,y,missingValue, /variables, filename=flnm_sav
endif ;sw_save



 print,' end pro restore_sza'
  end ;pro restore_sza
