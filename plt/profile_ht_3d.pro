;20110803! CAUTION! SH does not work!!!
;24 Mar 2011
;copied from plot_ptofile_ht.pro extracted plotting part only...
PRO  profile_ht_3d,plot_x,plot_y, title_hemi,mlat_title,ut_hr,lt_hr $
,plot_DIR,FLDIM_plot,mp_plot,sw_debug,sw_fort $
,sw_dif,sw_output2file,n_file


fac_window=0.50

n_ldct=39

y_min = 90.
y_max =220.;680.;  2.0E+04 ;18783.500 ;MAX(zn)
if ( sw_debug eq 1 ) then  print, "y_max=", y_max

if ( sw_fort eq 168L ) then begin
x_min=[-15.,   +100.,     0. ,     0.]
;x_max=[ 6.5,  +6500.,  +700. , +7000.]
x_max=[ 4.0,  +1300.,  +700. , +7000.]
endif else if ( sw_fort eq 167L ) then begin
x_min=[ 5.2,  + 60.,   -1. ,   0.]
x_max=[ 11.,  +180.,   +0.6, +14.]
endif


;get FLDIM and n_file
size_result=SIZE(plot_x)
if ( sw_debug eq 1 ) then  print,size_result
plot_type_max=size_result[1]
k_species    =size_result[2]
FLDIM        =size_result[3]
;n_file       =size_result[4] !this does not work when n_file=1
;if ( sw_debug eq 1 ) then  $
  print,plot_type_max,k_species," FLDIM=",FLDIM," n_file=",n_file


; set up plot parameters
i_plot=INTarr(n_file)
  i_plot[0]=0 ;file[0]
if ( n_file ge 2 ) then  i_plot[1]=1 ;file[1]
if ( n_file ge 3 ) then  i_plot[2]=2 ;file[2]
if ( n_file ge 4 ) then  i_plot[3]=3 ;file[3]
if ( n_file ge 5 ) then  i_plot[4]=4 ;file[3]
line_color=INTarr(n_file)
line_color[ i_plot[0]  ] =  55. ;blue
if ( n_file ge 2 ) then  line_color[ i_plot[1]  ] = 250. ;red
if ( n_file ge 3 ) then  line_color[ i_plot[2]  ] = 150. ;green
if ( n_file ge 4 ) then  line_color[ i_plot[3]  ] = 200. ;orange
if ( n_file ge 5 ) then  line_color[ i_plot[4]  ] = 100. ;orange

line_thick=[4.5,3.5,3.5,2.5,2.5]

;k_species = 4
;k_plot    =INTarr(k_species)
line_style=INTarr(k_species)

;title_var=['O+/H+/He+/Min+', 'Ti/Te/Tn', 'Vo+/Un/fo+', 'PHION/EHT']
if ( sw_dif eq 0 ) then begin
   if ( sw_fort eq 168L ) then  title_var=['O+/H+/He+/N+', 'Ti/Te', 'PHION/SUMION', 'EHT'] else $
   if ( sw_fort eq 167L ) then  title_var=['SL', 'GL/SZA', 'BM/GR', 'O/H+He/N2+O2/N4S']
endif else if ( sw_dif eq 1 ) then begin
   if ( sw_fort eq 168L ) then  title_var=['dif:O+/H+', 'dif:Ti/Te', 'PHION/SUMION', 'EHT'] else $
   if ( sw_fort eq 167L ) then  title_var=['SL', 'GL/SZA', 'BM/GR', 'O/H+He/N2+O2/N4S']
endif
;title_type='NTW'
line_style=[ 1,2,3,0 ]


!P.MULTI=[0,1,2,0,1]
;1:# of plot columns
;2:# of rows

;print, 'before',!P.BACKGROUND
!P.BACKGROUND=255
;print, 'after',!P.BACKGROUND

; plot height profile
  device, decomposed = 0 ,retain=2
;  window, 0 ,XSIZE=800,YSIZE=550
  window, 0 ,XSIZE=800*fac_window,YSIZE=750*fac_window
  loadct, n_ldct ;=39

  axis_color = 0.0 ;255.99 $
  char_size = 1.5

for plot_type=0,plot_type_max-1  do begin
;if ( plot_type eq 0 ) then   k_plot    =[ 0,1,3]  else $ ;O+  ;H+ ;Min+; He+ ;NNO
;if ( plot_type eq 1 ) then   k_plot    =[ 4,5,6]         ;Ti  ;Te ;Tn
;if ( plot_type eq 2 ) then   k_plot    =[ 4,5,6]         ;Vo+ ;Un ;o+flux
;if ( plot_type eq 3 ) then   k_plot    =[ 4,5,6]         ;phion ;eht


if ( sw_debug eq 1 ) then  print,'plot_type=',plot_type," x_min=", x_min[plot_type]," x_max=", x_max[plot_type]




i=0
title_plot=title_hemi+': mlat[deg]='+mlat_title+'  UT[hrs]='+STRTRIM( string(ut_hr[i_plot[i]], FORMAT='(f7.2)'), 1)+'  LT[hrs]='+STRTRIM( string(lt_hr[i_plot[i]], FORMAT='(f7.2)'), 1)

;plotting NH only
;if ( title_hemi eq 'NH' ) then begin
  j0=0L
  j1=FLDIM_plot[0]/2 -1L
;endif else if ( title_hemi eq 'SH' ) then begin
;  j0=FLDIM_plot[0]/2
;  j1=FLDIM_plot[0] -1L
;endif

;dbg print,'!debug plot_y',plot_y[j0:j1,i_plot[i]]
;frame only
if ( plot_type eq 0 ) then begin  ;densities
  plot, plot_x[plot_type,0,j0:j1,i_plot[i]], plot_y[j0:j1,i_plot[i]] $;, /xlog $ 
  ,xrange=[x_min[plot_type],x_max[plot_type]], xstyle=1  $
  ,yrange=[y_min,y_max], ystyle=1  $
  ,title=title_var[plot_type] $;+': '+title_plot  $
  ,linestyle = 0 $
  ,color=axis_color $
  ,charsize=char_size $
  ,/NODATA

endif else if ( plot_type ge 1 ) then begin  ;temperatures

  plot, plot_x[plot_type,0,j0:j1,i_plot[i]], plot_y[j0:j1,i_plot[i]] $
  ,xrange=[x_min[plot_type],x_max[plot_type]], xstyle=1  $
  ,yrange=[y_min,y_max], ystyle=1  $
  ,title=title_var[plot_type] $
  ,linestyle = 0 $
  ,color=axis_color $
  ,charsize=char_size $
  ,/NODATA

endif

for k=0,k_species-1 do begin
  for i=0,n_file-1 do begin
;d print,'check i=',i,n_file
;dbg20110803
;d   if i ge 3 then continue


;if ( title_hemi eq 'NH' ) then begin
  j0=0L
  j1=FLDIM_plot[i]/2-1L
;endif else if ( title_hemi eq 'SH' ) then begin
;  j0=FLDIM_plot[i]/2
;  j1=FLDIM_plot[i]-1L
;endif

if ( sw_debug eq 1 ) then  print, i, i_plot[i], line_color[i_plot[i]]
if ( sw_debug eq 1 ) then  print, plot_x[plot_type,k,j0:j1,i_plot[i]]
    oplot, plot_x[plot_type,k,j0:j1,i_plot[i]], plot_y[j0:j1,i_plot[i]] $ 
    ,linestyle = line_style[k] $
    ,color=line_color[i_plot[i]] $
    ,thick=line_thick[i]
 endfor ;i
endfor ;k

endfor ;plot_type=0,2 do begin

xyouts, 0.20, 0.5, title_plot  $
, charsize=1.0, charthick=1.0, color=axis_color, /norm, /noclip

if ( sw_output2file eq 'PNG' ) then begin
filename=plot_DIR+'profile_ht'+STRTRIM( string(sw_fort, FORMAT='(i3)'), 1)+'.mlat'+mlat_title+'_UT'+STRTRIM( string(ut_hr[0], FORMAT='(f7.2)'), 1)+'_'+title_hemi+'_mp'+STRTRIM( string((mp_plot+1), FORMAT='(i2)'), 1)+'dif.png'
output_png, filename
endif

;STOP


end ;profile_ht_3d
