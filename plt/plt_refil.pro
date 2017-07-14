;20120322UNDERCONSTRUCTION!!!
;date
;purpose: comparison with FLR mass density plots by Obana
;called from plt_ipe_sh.pro when plot_type=4
pro plt_refil   , mlat_deg, JMIN_IN,JMAX_IS, plot_z,mp_plot,n_read_max,ut_hr_save,fac_window, plot_UT,plot_UT_end,runID,plot_DIR,sw_debug

sw_save=1L
fac_den=1.0E-6
sw_refilling_rate=0L

;0number density
;sw_density=0
;1 mass density
sw_density=1L

if ( sw_density eq 0L ) then $
  y_title='Number Density [cm-3]' $
else if sw_density eq 1L then $
  y_title='MASS Density [amu /cc]'

;  DATA MASS/26.7616E-24,1.6726E-24,6.6904E-24,0.0/, [gram?]
;mass[4]=1.6726E-24 ;[gram]
;mass[5]=6.6904E-24 ;[gram]
;amu: A/16,1,4,0/ [amu]???
mass=fltarr(6)
mass[3]=16. ;[amu]
mass[4]=1.007 ;[amu]
mass[5]=4. ;[amu]



  FMT_l='(F7.2)'

jth_min=3
jth_max=6
den  =fltarr(jth_max+1,n_read_max)
ratio=fltarr(jth_max+1,n_read_max)
  for jth=jth_min,jth_max do begin          ;3:o+; 4:h+ ;5:he+; 6: total


;note20160907: h+, he+ filenames should not contain "+"
;US, 
     if mp_plot eq 2 then begin
        y_max=3000.             ;l=2.7
        locTitle='US'
;newZland 
     endif else if mp_plot eq 56 then begin
        y_max=5250.             ;l=2.74343
        locTitle='NZ'
     endif

        fac_ratio=1000.

          if jth eq 3 then $
        VarTitle='O' $
     else if jth eq 4 then $
        VarTitle='H' $
     else if jth eq 5 then $
        VarTitle='He'$
     else if jth eq 6 then $
        VarTitle='Total'


;l=3
;   if jth eq 4 then y_max=1200.  else $
;   if jth eq 5 then y_max= 110.; [cm-3] 

     print,'jth=',jth,VarTitle,'+'

imin=0
imax=0;2
for i=imin,imax do begin

if ( i eq 0 ) then begin
  lp=31 ;(1)US-canada
  ;lp=31  ;(2)NZ L=2.74343
  line_style=0
  line_thick=3.
endif else if ( i eq 1 ) then begin
  lp=25 ;L=3.94
  line_style=5
  line_thick=3.
endif else if ( i eq 2 ) then begin
  lp=23 ;L=5.13
  line_style=4
  line_thick=3.
endif

   if jth eq 3 then color_line=254 else $  ;o+
   if jth eq 4 then color_line=45 else $ ;H+
   if jth eq 5 then color_line=150 else $ ;He+
   if jth eq 6 then color_line=0  ;total
;color_line=50 ;50*(lp-lpstrt)/20
;d print, 'color_line', color_line

Re=6.3712E+03;km
r_ref=Re+90. ;km for reference ht for rcm
;assuming that latitude grid is identical for all LT sectors

midpoint = JMIN_IN[lp]-1 + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
theta = (90.- mlat_deg[  JMIN_IN[lp]-1 ])  *!PI / 180.  ;[deg]-->[rad]
sinthet = SIN( theta ) 
lval    = r_ref / ( Re * sinthet * sinthet )
print,'L-value=',lval,' mlat=',mlat_deg[jmin_in(lp)-1],' lp=',lp;,(jmin_in(lp)-1),(jmax_is(lp)-1)
;print,midpoint,mlat_deg[midpoint]


for n=0,n_read_max-1 do begin

   if sw_density eq 0 then begin
      den[jth.n]=plot_z[n,jth,0,midpoint]*fac_den ;m-3-->cm-3
      ratio[jth,n] =  plot_z[n,jth,0,midpoint] / plot_z[n,0,0,midpoint]
   endif else if sw_density eq 1 then begin

      if jth lt 6 then begin
         den[jth,n]=plot_z[n,jth,0,midpoint]*MASS[jth]*fac_den ;m-3-->cm-3

         ratio[jth,n] = ( plot_z[n,jth,0,midpoint]*MASS[jth] ) / ( $
                    plot_z[n,3,0,midpoint]*MASS[3] + $ ;o+
                    plot_z[n,4,0,midpoint]*MASS[4] + $ ;h+
                    plot_z[n,5,0,midpoint]*MASS[5]   ) ;he+
      endif else if jth eq 6 then begin                                 ;if jth lt 6 then begin

         den[jth,n]= $
            plot_z[n,3,0,midpoint]*MASS[3]*fac_den $
            + plot_z[n,4,0,midpoint]*MASS[4]*fac_den $
            + plot_z[n,5,0,midpoint]*MASS[5]*fac_den 
         
         ratio[jth,n]=-9999.9999
      endif

   endif                        ;sw_density eq 0 then begin
endfor                          ;n

if sw_refilling_rate eq 1L then calculate_refilling_rate, den, ut_hr_save, n_read_max,rate,sw_debug

if sw_debug eq 1 then  print,'Lval=',lval,' ratio=', ratio[jth,*]

x_min=0.;518400./3600.;86400./3600.;345600./3600.  ;plot_UT/3600.
x_max=24.;604800./3600.;950400./3600.;x_min+24.     ;plot_UT_end/3600. - plot_UT/3600.
y_min=0.







;dbg
print,jth,'check den: MAX', MAX(den[jth,*]),MIN(den[jth,*])
ut0=MIN(ut_hr_save)
print,'check ut MAX', MAX(ut_hr_save)
print,'check ut ut0=', ut0

if ( i eq imin AND jth eq jth_min ) then begin
  device, decomposed = 0 ,retain=2
  window, 0 ,XSIZE=800*fac_window,YSIZE=800*fac_window
  !P.BACKGROUND=255
   loadct, 0
   plot, (ut_hr_save-ut0), den[jth,*] $
; need to specify xmin,xmax
  ,xrange=[x_min,x_max], xstyle=1  $
  ,yrange=[y_min,y_max], ystyle=1  $
  ,title=runID+' '+' : L='+STRTRIM( string(lval, FORMAT=FMT_l), 1)+' '+locTitle $
  ,ytitle=y_title $
  ,linestyle = 0 $
  ,color=0. $ ;color_line 
  ,thick=line_thick $
, charsize=2.0, charthick=2.0 $
  ,/NODATA

endif ;if ( i eq imin AND jth eq jth_min ) then begin
   loadct, 39

   oplot, (ut_hr_save-ut0), den[jth,*] $
          ,linestyle =line_style $
          ,thick=line_thick $
          ,color=color_line 

;plot rati
   oplot, (ut_hr_save-ut0), (ratio[jth,*]*fac_ratio) $
          ,linestyle =2 $
          ,thick=line_thick $
          ,color=color_line 

endfor                          ;i=imin,imax

if ( jth eq jth_max ) then begin

   if mp_plot lt 10 then $
      FMT_mp='(i1)' $
   else if mp_plot lt 100 then $
      FMT_mp='(i2)'             ;$

   FILE_DISP=runID+VarTitle+'plus_L'+STRTRIM( string(lval, FORMAT=FMT_l), 1)+'_mp'+STRTRIM( string(mp_plot, FORMAT=FMT_mp), 1)+'_sw'+STRTRIM( string(sw_density, FORMAT='(i1)'), 1)+'.png'
   print,'L=',lval,' file_disp',file_disp
;  if ( sw_output2file eq 1 ) then  
   output_png, plot_DIR+FILE_DISP

   if sw_save eq 1L then begin
      FLNM_SAV=runID+VarTitle+'plus_L'+STRTRIM( string(lval, FORMAT=FMT_l), 1)+'_mp'+STRTRIM( string(mp_plot, FORMAT=FMT_mp), 1)+'_sw'+STRTRIM( string(sw_density, FORMAT='(i1)'), 1)+'.sav'
      print,'saving file=', flnm_sav
      save, ut_hr_save,ut0,den,ratio,lval,mp_plot,runID, /variables, filename=plot_DIR+FLNM_SAV
   endif                        ;sw_save eq 1
   
endif                           ; ( jth eq jth_max ) then begin 

endfor                          ;jth=4,5 do begin ;4:h+ ;5:he+




end ;pro plt_refil   , mlat_deg, JMIN_IN,JMAX_IS, plot_z
