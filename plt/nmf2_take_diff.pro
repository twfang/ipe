pro nmf2_take_diff

UserName='Raffaele.Montuoro'
RDIR00='/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/fig/'+UserName+'/'
RDIR0='124589'; 8x10;
RDIR1='128625'; 1x80;

;for utime=6,120,6 do begin
for utime=36,36 do begin
print,'utime=',utime

if utime lt 10 then $
  UTimeTitle=STRTRIM( string(utime, FORMAT='(i1)'),1 )+'.00' $
else if utime lt 100 then $
  UTimeTitle=STRTRIM( string(utime, FORMAT='(i2)'),1 )+'.00' $
else if utime lt 1000 then $
  UTimeTitle=STRTRIM( string(utime, FORMAT='(i3)'),1 )+'.00' 

;UTimeTitle='1.95'
print,'UTimeTitle=',UTimeTitle
jmax=0L;5L
VarTitle='NmF2';Un';On_m3'
htTitle='200.';350.'

;for j=0,jmax-1 do begin
;read continuous
flnm_sav0=RDIR00+'rt'+RDIR0+VarTitle+'UT'+UTimeTitle+'_ht'+htTitle+'.sav'
flnm_sav1=RDIR00+'rt'+RDIR1+VarTitle+'UT'+UTimeTitle+'_ht'+htTitle+'.sav'
print,'restoring file0=',flnm_sav0
   restore,flnm_sav0
;help
   print,'(0) utime=', ut_hr_disp


      var0=plot_zz
   
   max0=MAX(var0)
   min0=Min(var0)
   ave0=(max0+min0)*0.5
   print,'(0)nmf2:MAX=',max0,min0
   

;read restart
print,'restoring file1=',flnm_sav1
   restore,flnm_sav1
;help
   print,'(1) utime=', ut_hr_disp
   


      var1=plot_zz

   max1=MAX(var1)
   min1=Min(var1)
   print,'(1)nmf2:MAX=',MAX1,MIN1

   diff=var1-var0
   rat=diff/max0
   print,'diff=',max(diff,i),i,min(diff,j),j
   print,'rat=',max(diff)/max0, min(diff)/max0
   
;endfor                          ;j

;if max(diff) eq 0.0 then STOP
if max(diff) eq 0.0 then CONTINUE

   DEVICE, RETAIN=2, DECOMPOSED=0

fac_window=5.0
n_levels=100L
zmax=+0.04;1.
zmin=-0.04;-1.
x_min=-180.
x_min=+180.
Y_max= +80.0
Y_min= -Y_max

   WINDOW,0,XSIZE=1100*fac_window,YSIZE=1000*fac_window
;   loadct,39 ;n_ldct
   redblue
   contour,rat,plot_xx,plot_yy $
           ,/irregular $
           ,/fill $
           ,levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
           ,xrange=[X_min,X_max], /xstyle $
           ,yrange=[Y_min,Y_max], /ystyle $
           ,XTITLE=X_TITLE,YTITLE=Y_TITLE $
           ,TITLE='UT '+STRTRIM( string(ut_hr_disp, FORMAT='(F6.2)'),1 )+' MIN='+STRTRIM(STRING( MIN(rat), FORMAT='(E11.3)'),1)+' MAX='+STRTRIM(STRING( MAX(rat), FORMAT='(E11.3)'),1) $
           ;,POSITION=[X0,Y0,X1,Y1] $
           ,COLOR=text_color $
           ,charsize=char_size,charthick=char_thick $
           ,MAX_VALUE= MAX_xymin ;$
           ;,/OVERPLOT

   charsize_colorbar=4.0
   format_colorbar='(f7.2)'
   position=[0.30, 0.08, 0.8, 0.09] ;for horizontal bar
   color=0
   COLORBAR, BOTTOM=bottom, CHARSIZE=charsize_colorbar, COLOR=color, DIVISIONS=divisions $
             , FORMAT=format_colorbar, POSITION=position, MAXRANGE=zmax,MINRANGE=zmin $
             , NCOLORS=ncolors,TITLE=title,VERTICAL=vertical,TOP=top,RIGHT=right $
             , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
             , _EXTRA=extra, INVERTCOLORS=invertcolors, TICKNAMES=ticknames
Filename_png=RDIR00+'diff_'+RDIR0+'_'+RDIR1+'UT'+UTimeTitle+'.png'
print,Filename_png
output_png,Filename_png

endfor ;utime=6,120,6 do begin

print, "nmf2_take_diff finished!"
end;pro nmf2_take_diff
