;20141007 added interpolation for smoother plot and TEC/nmf2/hmf2 calculations 
;TODO
;include parallel V for debug purpose!!!
;20110203: copied from /home/naomi/sandbox/SUPIM/plot_idl/contour_result.pro

 PRO  ctr_lat_ht_tec $
, in2d,is2d,z_km,mlat_deg  $             ;input 
, plot_z,plot_VEXB,n_read   $ ;input
, uthr, plot_DIR, title_res $
;,rundate,title_test $
,sw_debug, title_hemi,sw_anim,mpstart,mpstop,mpstep, lt_hr, glon_deg, fac_window $
, sw_output2file, TEST, TEST2, TEST1 $
, VarType_min $;=8L
, VarType_max $;=8L ;PAR-1
, VarType_step $;=4L
, luntmpN ;
;
;VarType_min=3L
;VarType_max=3L ;PAR-1
;VarType_step=4L

print,' mpstart', mpstart,'mpstop',mpstop



sw_plot_tec=1
if ( sw_plot_tec eq 1 ) then begin
   sw_output_nmf2noon=0
   nct_max=44000L
   x_save=fltarr(nct_max)
   y_save=fltarr(nct_max)
   z_save=fltarr(nct_max)
;dbg20160127 this switch is used to output .txt faster without plotting.    
   sw_SaveTec=1
   if ( sw_SaveTec eq 1 ) then $
      sw_NoPlot=1 $
   else $
      sw_NoPlot=0
endif

fac_ne = 1.0E+12
sw_arw_vpara=0
fac_arw_para=2.0

n_read0=0L
sw_plot_grid=0L  ;1: plot grid only
sw_arrow_exb=0L
reference_arrow=40L  ;m/s
;factor_arrow=5.
factor_arrow=5.

lp_step_arrow=7
;lpmax_perp_trans=37;151-1
lpmax_perp_trans=149
;mpstart=mp_plot-5;0
;mpstop=mp_plot;0
;mpstep=1
;debug
for i=mpstart, mpstop  do print,' mp', (i+1),' LT',lt_hr[i]

HTmin=90.  ;min(yy)   ;75.   ;400. ;
HTmax=1040.;700.;190. ;1.000000E+03;700.; 
; plot range
if ( title_hemi eq 'NH' ) then begin
  gLATmax=+90.;+90.;-10.;
  gLATmin=+60.;+50.;-gLATmax;-27.; 
endif else if ( title_hemi eq 'SH' ) then begin
  gLATmax=-48.;+90.;-10.;
  gLATmin=-78.;+50.;-gLATmax;-27.; 
endif else if ( title_hemi eq 'glb' ) then begin
  gLATmax=+60.;-10.;
  gLATmin=-gLATmax;-27.; 
endif else if ( title_hemi eq 'eq' ) then begin
  gLATmax=-15.;+90.;-10.;
  gLATmin=-33.;-gLATmax;-27.;
;  HTmax=1.001E+03 
endif




;title_test='trans'
sw_dif=0L
device_type='png' ;ps';'


size_result = SIZE(in2d)
NLP=size_result[1]
if ( sw_debug eq 1 ) then  print, 'NLP=', NLP
size_result = SIZE(z_km)
NPTS2D=size_result[1]
if ( sw_debug eq 1 ) then  print, 'NPTS2D=', NPTS2D
size_result = SIZE(plot_z)
NPAR = size_result[2]




N_LDCT=39;33
;lp_strt=1;28-1;  0+1 ;58;0;63 ;1-1L
lp_strt=0+1 ;58;0;63 ;1-1L
;lp_stop=NLP-1-1 ;138L;
lp_stop=NLP-2;138L;



;if ( uthr le 0.00 ) then $
;  time_string='00'+STRTRIM(STRING( uthr-24.0, FORMAT='(F8.4)'),1)+'UT' $
;else 
if ( uthr lt 1.0E+01 ) then $
  time_string='00'+STRTRIM(STRING( uthr, FORMAT='(F8.4)'),1)+'UT'  $
else if ( uthr lt 1.0E+02 ) then $
  time_string='0'+STRTRIM(STRING( uthr, FORMAT='(F8.4)'),1)+'UT'  $
else $ ;if ( uthr ge 1.0E+02 ) then $
  time_string=STRTRIM(STRING( uthr, FORMAT='(F8.4)'),1)+'UT'

for mp=mpstart,mpstop, mpstep do begin
if ( sw_debug eq 1 ) then  print, 'mp=',mp

;042004:
;if ( elhr le 0.00 ) then $
FileID=time_string+'_'+STRTRIM(STRING( lt_hr[mp], FORMAT='(F6.2)'),1)+'LT'+'_mp'+STRTRIM(STRING( (mp+1), FORMAT='(i3)'),1)+'.'+TEST ;mp+1:ipe convention

VarTitle=[ $
'Ne',$
'Te',$ 
'Ti',$ 
'O+','H+','He+','N+','NO+','O2+','N2+','O+2D','O+2P' $
,'Vo+' $
;'o+flux'$
;'hr4',$
]

VarUnit=[ $
'[log!D10!N cm-3]',$
'[K]', $ ;Te
'[K]', $ ;Ti
'[log!D10!N cm-3]','[log!D10!N cm-3]','[log!D10!N cm-3]','[log!D10!N cm-3]','[log!D10!N cm-3]','[log!D10!N cm-3]','[log!D10!N cm-3]','[log!D10!N cm-3]','[log!D10!N cm-3]' $
;'[m/s]',$
;,'[cm2 s-1]' $
;'[J/kg/s]', $;hrate
,'[m/s]'$
]

X_Title='MAGNETIC LATITUDE [deg.]'
Y_Title='ALTITUDE [km]'

;RE=6.3712E+03  ;[km]


lp=0L
ihem=0L   ;0:north; 1:south
istrt=0L
istop=-1L
istep=0L
X=dblarr(4)
Y=dblarr(4)

if ( sw_dif eq 0 ) then begin

   ARY_min0=[ $
0.,$  ; X1.0E+12
;9.996,$   ;alog10
170.,$
170., $
0., $ ;o+
0., $ ;h+
0., $ ;he+
0., $ ;n+
0., $ ;no+
0., $ ;o2+
0., $ ;n2+
0., $ ;o+2D
0., $ ;o+2P
          -10.  $   
         ] 

ARY_max0=[ $
;6.2, $;
2.377,$ ;7.,$ ;X1.0E+12
;12.46,$ ;7.,$ ;log
;4.,$
;800. ,$
1400. ,$
1400. ,$
1.28, $ ;o+
1.9E-3, $ ;h+
5.7E-4, $ ;he+
8.7E-3, $ ;n+
2.5E-2, $ ;no+
5.0E-3, $ ;o2+
1.5E-3, $ ;n2+
1.9E-3, $ ;o+2D
6.8E-5, $ ;o+2P
          +10.     $  ;vel [m s-1]
        ]
endif else if ( sw_dif eq 1 ) then begin
ARY_min0=[ $
           3., 3., 3., 3. $  ;densities
        ,  -1.940E+2        ,  -1.940E+2  $      ;To+;Te
        ,  -1.0E+13  $;flux [cm2 s-1]
        ,  -3500.  $   
         ] 

ARY_max0=[ $
           7., 7., 7., 7.  $ ;densities
        ,  +1.253E+03,  +1.253E+03  $    ;To+,Te
        ,  +1.0E+13  $
        ,  +3500.     $  ;vel [m s-1]
        ]
endif ;( sw_dif eq 1 ) then begin

col_min=0.   +0.3
col_max=256. -0.3
N_LVLs=100  ;intarr(2)       ; number of contour levels
;*** Load color table
LOADCT, N_LDCT



; window size
X_SIZE=27. ;9.7
Y_SIZE=21.
X0= 3.0 ;[cm]
Y0= 3.0 ;
dX= 23.0
dY= 15.0

;color bar
dX1=  7.5
dY1=  0.2
X1=X0+dX-dx1 ;+1.5    ;7.5   ;8.   ;18.1
Y1=Y0+dY+1.3 ;-0.5     ;+ 1.6   ;1.3   ;18.8

xx=dblarr(npts2D)
yy=dblarr(npts2D)
;; plot GLAT v.s. ALT Flux-Tube Distribution
for lp=0,NLP-1 do begin   ;ifl-->lp

  in = IN2D[lp]-1L
  is = IS2D[lp]-1L 
;  print, 'mp',mp,' lp',lp,' in',in,' is',is

  istrt=istop +1
  istop=istrt+(is-in) ;istrt +(iseb(ifl)-ineb(ifl))

;  print,'istrt=',istrt,' istop=',istop

; set up values to be plotted
;!!!!UNDERCONSTRUCTION!!!
  xx[istrt:istop] = mlat_deg[in:is]  ;glatd(ineb(ifl):iseb(ifl), ifl)
  yy[istrt:istop] = z_km[in:is]  ;gpz(  ineb(ifl):iseb(ifl), ifl)
endfor  ;lp=lp_strt,lp_stop do begin



for VarType=VarType_min , VarType_max,  VarType_step   do begin
if ( sw_debug eq 1 ) then  print,'plotting ',VarTitle(VarType)

if ( sw_plot_tec eq 1 ) then begin
nct=-1L
x_save[*]=0.0
y_save[*]=0.0
z_save[*]=0.0
endif

MainTitle=VarTitle(VarType)+' '+VarUnit(VarType)
FILE_DISP=plot_DIR+device_type+'/'+title_hemi+'/'+FileID+'_'+VarTitle(VarType)+'_'+title_res+'.'+title_hemi+'.'+device_type
if ( sw_debug eq 1 ) then  print, file_disp
if ( device_type eq 'ps' ) then begin



  SET_PLOT, 'ps' ;'x'
   ; Open a file in the current directory 
   ; to contain the PostScript plot.
  DEVICE, BITS_PER_PIXEL=24,/COLOR, $
  FILENAME=FILE_DISP, /LANDSCAPE,    $
  SET_FONT='Courier',       /TT_FONT, $
  XSIZE=X_SIZE, YSIZE=Y_SIZE, XOFFSET=0.0, YOFFSET=X_SIZE
endif else  if ( device_type eq 'png' ) then begin

;  !P.MULTI=[4,4,2,0,1]
;  ;1:# of plot columns
;  ;2:# of rows

if ( sw_debug eq 1 ) then  print, 'before',!P.BACKGROUND
  !P.BACKGROUND=0 ;255
if ( sw_debug eq 1 ) then  print, 'after',!P.BACKGROUND

; plot height profile
if ( sw_NoPlot eq 0 ) then  device, decomposed = 0 ,retain=2
if ( sw_NoPlot eq 0 ) then  window, 0 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

;t  axis_color = 0.0 ;255.99 $


endif

LOADCT, N_LDCT

;if ( sw_plot_grid eq 1 ) then begin ;20120328
;LOADCT, 0
if ( sw_NoPlot eq 0 ) then Plot, xx(0:istop), yy(0:istop)     $
, Xstyle = 1, Xrange = [ gLATmin, gLATmax]  $
, Ystyle = 1, Yrange = [ HTmin, HTmax]      $
, TITLE = MainTitle+'   '+FileID, SUBTitle =' ' $   ;FileID $
, XTITLE = X_Title,    YTITLE = Y_Title  $
, PSYM =  3,  SYMSIZE=1.0  $
;, Color = col_min  $
;, CharSize = 1.5 $
;, THICK    = 1.0 $
, Pos = [X0/X_SIZE, Y0/Y_SIZE, (X0+dX)/X_SIZE, (Y0+dY)/Y_SIZE]  $
,/NODATA

;if ( sw_output2file eq 1 ) then begin
;  print,'output to file=',plot_DIR+'ipe_grid.'+device_type
;  output_png,plot_DIR+'ipe_grid.v2.'+device_type
;endif
;LOADCT, N_LDCT
;RETURN

;endif ;( sw_plot_grid ne 1 ) then begin ;20120328


;tmp20110414: plot grid only
;if ( sw_plot_grid ne 1 ) then begin

ARY_minZ=ARY_max0(VarType)
ARY_maxZ=ARY_min0(VarType)
;print, '0ARY_minZ', ARY_minZ
;print, '0ARY_maxZ', ARY_maxZ

for lp=lp_strt , lp_stop do begin

;d if ( sw_debug ) then $
;if ( (lp+1) ge 148 ) AND ( (lp+1) le 155 ) then $
;d  print,'lp=', lp, in2d[lp], mlat_deg[ in2d[lp] ]

; midpoint = JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
  midpoint =      IN2D[lp]+ (     IS2D[lp] -    IN2D[lp]    )/2  -1 

;d if ( sw_debug ) then $
;if ( (lp+1) ge 148 ) AND ( (lp+1) le 155 ) then $
;d print,(lp+1),midpoint,'mlat',mlat_deg[ midpoint ]  ,'midpoint z', z_km[ midpoint ],' max z',MAX( z_km[ IN2D[lp]:IS2D[lp] ] )


  for ihem=0,1   do begin

    istrt=midpoint
    if ( ihem eq 0 ) then begin ;North
      istep= -1
      istop=in2d[lp]-1 ;-istep 
      factor=+1.0
    endif else if ( ihem eq 1 ) then begin   ;South
      istep= +1
      istop=is2d[lp]-1 ;-istep 
      factor=-1.0
    endif


;if ( mlat_deg[ is2d[lp]-1 ] ge -66. ) AND ( mlat_deg[ is2d[lp]-1 ] lt -62. ) THEN  begin
;print,'check mlat!', is2d[lp], lp ,  mlat_deg[ is2d[lp]-1 ]
;stop
;endif


;if ( lp ge 50 ) then $
  ;dXX=0.7; (mlat_deg[ in2d[lp-1] ]-mlat_deg[ in2d[lp] ])*0.95*factor ;nm20121130
  dXX=(mlat_deg[ in2d[lp-1] ]-mlat_deg[ in2d[lp] ])*0.9*factor ;nm20121130
;  dXX=(mlat_deg[ in2d[lp-1] ]-mlat_deg[ in2d[lp] ])*0.95*factor
;else $
;  dXX=1.00

    for ipts=istrt,istop,istep   do begin
;dYY=(z_km[ipts+1]-z_km[ipts] )*1.0 ;
;dYY=5.;1.6;(z_km[ipts+1]-z_km[ipts] )*0.9 ;nm20121130
dYY=(z_km[ipts+1]-z_km[ipts] )*0.9 ;nm20121130

if ( z_km(ipts  ) gt HTmin ) and ( z_km(ipts) lt (HTmax-dYY) ) then begin
if ( mlat_deg(ipts ) gt gLATmin ) and ( mlat_deg(ipts) lt gLATmax ) then begin

;dbg20140828
;if ( mlat_deg[ipts] gt -24. ) AND  ( mlat_deg[ipts] lt -22. )  AND (z_km[ipts] gt 340. ) AND (z_km[ipts] lt 344. ) then $
;print,'lp', lp,' mlat',mlat_deg[ipts],' z_km',z_km[ipts],plot_z[n_read,VarType, 0,ipts]*1.0E-12

Xa=mlat_deg(ipts) ;-dXX*.5
Xb=mlat_deg(ipts) ;-dXX*.5   ;glatd(ipts  ,ifl)+dLAT  ;
Xc=mlat_deg(ipts)+dXX ;*.5  ;Xb  ;
Xd=mlat_deg(ipts)+dXX ;*.5  ;Xa  ;
Ya=z_km(ipts) ;-dYY*.5
Yb=z_km(ipts)+dYY ;*.5 ;Ya ;
Yc=z_km(ipts)+dYY ;*.5 ;Ya+dHT ;
Yd=z_km(ipts) ;-dYY*.5   ;Yc     ;

;nm20141006
if ( sw_plot_tec eq 1 ) then begin
nct=nct+1
x_save[nct]= mlat_deg[ipts]
y_save[nct]= z_km[ipts]
z_save[nct]= ALOG10( plot_z[n_read,VarType, 0,ipts] ) ;[m-3]
endif

 X=[Xa, Xb, Xc, Xd]  ;glat [deg]
 Y=[Ya, Yb, Yc, Yd]  ;altitude [km]

if ( VarType eq 0 ) OR ( VarType ge 3 ) then begin

  density = $
    plot_z[n_read,VarType, 0,ipts] ;t * 1.0E-6  ;m-3 --> cm-3

  if ( density gt 0.0 ) then $
     Value= ALOG10( density ) $
  else $ 
     Value= ALOG10( 0.1 )
;20131204
Value=plot_z[n_read,VarType, 0,ipts] / fac_ne

endif else if ( VarType eq 1 ) or ( VarType eq 2 ) then $ ;Te/i
  Value = plot_z[n_read,VarType,0,ipts] ;



; save actual MIN & MAX values
if ( Value lt ARY_minZ ) then  ARY_minZ=Value

if ( Value gt ARY_maxZ ) then begin
 ARY_maxZ=Value
; print,'ARY_maxZ',ARY_maxZ,' z_km',z_km(ipts),' mlat',mlat_deg(ipts)
endif

ColorData=(col_max-col_min)*(Value-ARY_min0(VarType))/(ARY_max0(VarType)-ARY_min0(VarType)) + col_min


;061904: for fancy looking
if ( ColorData lt col_min ) then  ColorData=col_min
if ( ColorData gt col_max ) then  ColorData=col_max

;091604: debug
;091604: if ( (ipts MOD 2)  eq 0 ) then $
if ( sw_NoPlot eq 0 ) then $
POLYFILL, X , Y       $      ;[, Z]]
, FILL_PATTERN=0      $      ; solid fill
, COLOR=ColorData    ;$
;[, IMAGE_COORD=array]
;[, /IMAGE_INTERP] 
;[, /LINE_FILL] 
;[, PATTERN=array] 
;[, SPACING=centimeters] 
;[, TRANSPARENT=value]
;[, CLIP=[X 0 , Y 0 , X 1 , Y 1 ]] 
;[, /DATA | , /DEVICE | , /NORMAL] 
;[, LINESTYLE={0 | 1 | 2 | 3 | 4 | 5}] 
;[, /NOCLIP] 
;[, ORIENTATION=ccw_degrees_from_horiz] 
;[, /T3D] 
;[, THICK=value] 
;[, Z=value]

if ( sw_arw_vpara eq 1 ) AND $
 ( (lp MOD lp_step_arrow) eq 2 ) AND $
 ( (ipts MOD 10) eq 0 ) then begin
x2= (mlat_deg[ipts+1]-mlat_deg[ipts])*(mlat_deg[ipts+1]-mlat_deg[ipts])
y2= (z_km[ipts+1]-z_km[ipts])*(z_km[ipts+1]-z_km[ipts])
dr=SQRT( x2 + y2 )
cosI=SQRT( x2 ) / dr
sinI=SQRT( y2 ) / dr
;
X0_arrow = mlat_deg[ipts]
Y0_arrow = z_km[ipts]
Vpara    = plot_z[n_read,VarType,mp,ipts]*fac_arw_para
dX_arrow = Vpara *cosI
dY_arrow = Vpara *sinI
X1_arrow = X0_arrow + dX_arrow
Y1_arrow = Y0_arrow + dY_arrow

LOADCT, 0
head_thick=2.0
head_size=!D.X_SIZE / 64.*0.5 
if ( sw_debug ) then $
 print, 'head_size', head_size
;
  ARROW, X0_arrow, Y0_arrow, X1_arrow, Y1_arrow  $
, /DATA       $ ;???
, /NORMALIZED $
, HSIZE=head_size $
, COLOR=col_max  $ ;index] $
, HTHICK=head_thick $ ;value] $
;, /SOLID] $
, THICK=3 ;value]
LOADCT, N_LDCT
endif ;( sw_arw_vpara ) then begin


endif ;( glatd(ipts  ,ifl-1) gt gLATmin ) and ( glatd(ipts  ,ifl-1) lt gLATmax) thenbegin
endif  ;( gpz(ipts  ,ifl-1) gt HTmin ) then begin
endfor  ;ipts=istrt,istop,istep   do begin
endfor  ;ihem=1,1   do begin



;draw arrow of EXB drift at midpoint

if ( z_km[midpoint] le HTmax ) AND ( z_km[midpoint] ge HTmin ) then  begin

if ( sw_arrow_exb eq 1 ) AND (lp le lpmax_perp_trans) then begin
if ( lp le 46 ) $ ;mlat-30 boundary
OR ( lp eq lpmax_perp_trans ) $ ;transport boundary
OR ( (lp MOD lp_step_arrow) eq 2  ) then begin

;print, 'check exb!',n_read,mp,lp,plot_VEXB[n_read,mp,lp], z_km[midpoint]

X0_arrow = mlat_deg[midpoint]
Y0_arrow = z_km[midpoint]
dX_arrow = 0.
dY_arrow = plot_VEXB[n_read,mp,lp] * factor_arrow
X1_arrow = X0_arrow + dX_arrow
Y1_arrow = Y0_arrow + dY_arrow

LOADCT, 0
head_thick=2.0
head_size=!D.X_SIZE / 64.*0.5 
if ( sw_debug ) then $
 print, 'head_size', head_size
  ARROW, X0_arrow, Y0_arrow, X1_arrow, Y1_arrow  $
, /DATA       $ ;???
, /NORMALIZED $
, HSIZE=head_size $
, COLOR=col_max  $ ;index] $
, HTHICK=head_thick $ ;value] $
;, /SOLID] $
, THICK=3 ;value]
LOADCT, N_LDCT
endif ;( (lp MOD lp_step_arrow) eq 0  ) then begin
endif ;( sw_arrow_exb ) then begin

endif ;( z_km[midpoint] lt ) then  begin

endfor  ;lp=lp_strt,lp_stop do begin


;draw arrow of reference for 50m/s
if ( sw_arrow_exb eq 1 ) then begin
X0_arrow = gLATmax +  2.0
Y0_arrow = HTmin   - 60.
dX_arrow = 0.
dY_arrow = reference_arrow * factor_arrow
;dY_arrow = reference_arrow * fac_arw_para
X1_arrow = X0_arrow + dX_arrow
Y1_arrow = Y0_arrow + dY_arrow

LOADCT, 0
head_thick=2.0
head_size=!D.X_SIZE / 64.*0.5 
if ( sw_debug ) then $
 print, 'head_size', head_size
  ARROW, X0_arrow, Y0_arrow, X1_arrow, Y1_arrow  $
, /DATA       $ ;???
, /NORMALIZED $
, HSIZE=head_size $
, COLOR=col_max  $ ;index] $
, HTHICK=head_thick $ ;value] $
;, /SOLID] $
, THICK=3 ;value]

xyouts, (X0+dX-0.4)/X_SIZE, (Y0-0.7)/Y_SIZE, STRTRIM(STRING(reference_arrow , FORMAT='(I2)'),1)+'m/s'  $
, charsize=1.0, charthick=1.0, /norm, /noclip
LOADCT, N_LDCT

endif ;( sw_arrow_exb ) then begin


text_color=255.;depends on the background
if ( sw_NoPlot eq 0 ) then $
Draw_Colorbar, ARY_min0(VarType), ARY_max0(VarType), N_LVLs $
, col_min, col_max, X1, Y1, dX1, dY1, X_SIZE, Y_SIZE, VarType, text_color

;LT
;xyouts, (X0-2.)/X_SIZE, (Y0+dY+2.45)/Y_SIZE, STRTRIM(STRING( zthr, FORMAT='(F5.2)'),1)+'LT'  $
;, charsize=1.0, charthick=1.0, /norm, /noclip


print,'MIN',ary_minz,' MAX',ary_maxz
; add MIN & MAX values
if ( sw_NoPlot eq 0 ) then $ 
xyouts, (X0-0.8)/X_SIZE, (Y0+dY+0.6)/Y_SIZE, 'MIN='+STRTRIM(STRING( ARY_minZ, FORMAT='(E11.3)'),1)+' MAX='+STRTRIM(STRING( ARY_maxZ, FORMAT='(E11.3)'),1)  $
, charsize=1.0, charthick=1.0, /norm, /noclip


;STOP

if ( sw_debug eq 1 ) then  print, VarTitle(VarType), '  ARY_minZ=', ARY_minZ, ' ARY_maxZ=', ARY_maxZ
;endif ;( sw_plot_grid ne 1 ) then begin


    ; Close the "OUTPUT_DEVICE".
if ( device_type eq 'ps' ) then $
  DEVICE, /CLOSE $
else if ( device_type eq 'png' ) then begin

  if ( sw_anim eq 1 ) then FILE_DISP=plot_DIR+'anim/'+STRTRIM( string( (n_read+1), FORMAT='(i3)'), 1)+'.'+device_type

  if ( sw_output2file eq 1 ) then  output_png, FILE_DISP
endif


if ( sw_plot_tec eq 1 ) then begin
x3=x_save[0:nct]
y3=y_save[0:nct]
z3=z_save[0:nct]

; triangulate prosedure                                                         
TRIANGULATE,X3,Y3,Tr
if ( sw_debug eq 1 ) then help, Tr

d_x=1.0 ;deg                                                                         
d_y=20.;km

trigrid_result=trigrid(x3,y3,z3,tr  $
, [d_x,d_y], [gLATmin,HTmin,gLATmax,HTmax] $
;, NX=12, NY=24  $                                                              
)
size_result=size(trigrid_result)
if ( sw_debug eq 1 ) then  print, 'check trigrid_result', size_result , min(trigrid_result),max(trigrid_result)

nx=size_result(1) ;121
ny=size_result(2) ;49

x_dsp0=findgen(nx)*d_x + gLATmin
;for i=0,nx do begin
;print,' x_dsp0', i,x_dsp0[i]
;endfor
;stop

y_dsp0=findgen(ny)*d_y +   HTmin
;print,' y_dsp0', y_dsp0
z_dsp0=(10^(trigrid_result)) /fac_ne
;z_dsp0=trigrid_result ;alog10

zmax=MAX(z_dsp0)
zmin=MIN(z_dsp0)
print, 'Ne^10-12 zmin=',zmin,' zmax=',zmax
zmax=ARY_max0[VarType]
zmin=ARY_min0[VarType]
if ( sw_NoPlot eq 0 ) then   device, decomposed = 0 ,retain=2
if ( sw_NoPlot eq 0 ) then  window, 1 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

X_max=+55.
X_min=-X_max
Y_min=90.
Y_max=1000.
char_size = 1.8
char_thick = 2.
line_thick = 3.
text_color=0.

 !P.BACKGROUND=255

if ( sw_NoPlot eq 0 ) then contour,z_dsp0,x_dsp0,y_dsp0 $
,/fill $
, levels=findgen(N_LVLs)*(zmax-zmin)/float(N_LVLs-1) +zmin $
, xrange=[X_min,X_max], xstyle=1  $
, yrange=[Y_min,Y_max], ystyle=1  $
, XTITLE = X_Title $;'mlat [deg]' $
, YTITLE = Y_Title $;'ht [km]' $ 
,TITLE = 'Electron Density [10!U12!N  m!U-3!N]' $
, Pos     =[X0/X_SIZE, Y0/Y_SIZE, (X0+dX)/X_SIZE, (Y0+dY)/Y_SIZE]  $
, charsize = char_size, charthick = char_thick $
, COLOR=text_color

X1 = X1 - 0.3
Y1 = Y1 + 0.5
if ( sw_NoPlot eq 0 ) then Draw_Colorbar, zmin, zmax, N_LVLs $
, col_min, col_max, X1, Y1, dX1, dY1, X_SIZE, Y_SIZE, VarType, text_color

if ( sw_NoPlot eq 0 ) AND ( sw_output2file eq 1 ) then output_png, 'fig2_'+FileID+'_'+TEST2+'_'+TEST1+'.png'

TEC=fltarr(nx)
nmf2=fltarr(nx)
hmf2=fltarr(nx)
TEC(*)=0.0e0
for ix=0,nx-1  do begin  ;glat                                                  
  for iy=0,ny-1  do begin  ;htkm                                                

  ; log10 -> Ne[m-3]                                                                    
   ne1=10^(trigrid_result(ix,iy))

   if ( ne1 gt nmf2[ix] ) then begin
     nmf2[ix]=ne1
     hmf2[ix]=y_dsp0[iy]
   endif
   TEC[ix] = TEC[ix] + ne1*(d_y*1.0e+3)  ;[m-3 * m] =[m-2]                        
  endfor ;iy=0,nx-1  do begin  ;htkm                                            
;d print,'check TEC', ix, x_dsp0[ix],tec[ix]
endfor ;ix=0,nx-1  do begin 

TECU=1.0E16 ;[m-2]



;tmp20160116 !P.BACKGROUND=0
if ( sw_NoPlot eq 0 ) then  device, decomposed = 0 ,retain=2
if ( sw_NoPlot eq 0 ) then  window, 2 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

tec0=tec/TECU
;print,'MAX TEC',MIN(tec0),MAX(tec0)
;print,"tec0=[$"
;print, tec0, format='(5(",",f12.7),"$")'
;print, "]"
if ( sw_NoPlot eq 0 ) then Plot, x_dsp0, tec0 $
;, yrange=[ 0., 120. ], ystyle=1  $ ;20141119
;, yrange=[ 0., 65. ], ystyle=1  $
, yrange=[ 0., 80. ], ystyle=1  $
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50. , +50. ], xstyle=1  $
, TITLE='TEC: '+FileID $
,color=text_color, THICK =line_thick, charsize=char_size, charthick=char_thick

;v4 windX1.5
 tec4=fltarr(nx)
; tec4 =[$
;  0.0000000,   9.6984081,  10.1598568,  10.5176029,  11.0120668$
;,  11.6842537,  12.4723835,  13.3431339,  14.1974840,  15.4216709$
;,  16.9979343,  18.9607029,  21.1281395,  23.5545616,  26.1064434$
;,  28.9344692,  32.2859535,  35.7561836,  39.4601250,  43.3099480$
;,  46.9687157,  50.3004646,  53.5319557,  56.6343002,  59.3193855$
;,  61.6547623,  63.8377724,  65.4999619,  66.6006851,  67.6338882$
;,  68.2777863,  67.9558334,  67.7838745,  65.9939880,  59.4125977$
;,  51.0568848,  45.1089287,  41.5212173,  39.1856842,  37.5787277$
;,  36.5636482,  35.9289970,  35.5760078,  35.4573479,  35.6043968$
;,  35.9724464,  36.4794807,  37.0107155,  37.4963799,  37.8409615$
;,  38.0100861,  37.9833641,  37.8214645,  37.6032562,  37.3945465$
;,  37.2615013,  36.9876175,  36.8452797,  36.7807350,  36.7037201$
;,  36.5725822,  37.0751801,  38.1122932,  39.3375511,  40.8677864$
;,  42.5527878,  44.0995178,  45.6302147,  46.9605560,  47.9323196$
;,  48.2332764,  47.7684593,  46.4083900,  44.4957581,  42.4437027$
;,  40.2722321,  38.5130653,  36.9459991,  35.4471970,  34.2024918$
;,  33.1776009,  32.5048447,  32.0645981,  31.9518661,  32.0242958$
;,  32.2455597,  31.2924824,  26.5000000,  19.5464020,  16.3490219$
;,  14.5790977,  13.6063013,  13.0062675,  12.5200682,  12.1970072$
;,  11.9498720,  11.6956863,  11.4751863,  11.3143749,  11.1666451$
;,  11.0567837,  10.9933405,  10.9201870,  10.8451109,  10.8060074$
;,  10.7830391,  10.7722282,  10.7738829,  10.7823267,  10.7712002$
;,  10.7501488,  10.7175722,  10.7015343,  10.7311926,  10.6826897$
;,  10.6114788,  10.5184202,  10.4382181,  10.4025517,  10.3749847$
;,   0.0000000$
; ]
nmf24=fltarr(nx)
hmf24=fltarr(nx)
tec5=fltarr(nx)
nmf25=fltarr(nx)
hmf25=fltarr(nx)
tec6=fltarr(nx)
nmf26=fltarr(nx)
hmf26=fltarr(nx)
;dbg20160126 tests for GRL reply
;get_v4 $
;get_v20,time_string $
;get_v7,time_string $
;,tec4,nmf24,hmf24 $
;,tec5,nmf25,hmf25 $
;,tec6,nmf26,hmf26

if ( sw_NoPlot eq 0 ) then oPlot, x_dsp0, tec4, linestyle=5 ,color=col_min, THICK =line_thick
;
;v5 windX2
; tec5=fltarr(nx)
; tec5 =[$
;   0.0000000,   8.5055103,   8.8729372,   9.1701126,   9.5792685$
;,  10.1371927,  10.8040371,  11.5533123,  12.3200064,  13.4604206$
;,  14.9593668,  16.9990807,  19.4007511,  22.2253170,  25.3013592$
;,  28.7844200,  32.9525909,  37.3187675,  41.9807968,  46.7886620$
;,  51.3016624,  55.3511696,  59.1525345,  62.6950035,  65.6605225$
;,  68.0910568,  70.2884903,  71.8594513,  72.6758575,  73.2440033$
;,  73.1018448,  70.9367676,  66.7925720,  58.6078491,  47.7719040$
;,  40.5865631,  36.7204247,  34.5876884,  33.2943268,  32.5857201$
;,  32.3186302,  32.1868629,  32.1772041,  32.3115883,  32.6651802$
;,  33.1723061,  33.7973404,  34.4066162,  34.9293365,  35.2924309$
;,  35.4707222,  35.4450836,  35.2884064,  35.0857391,  34.9029007$
;,  34.8007622,  34.5829811,  34.5121078,  34.5400887,  34.5700836$
;,  34.4571304,  34.9368553,  35.9108200,  37.0806656,  38.5330162$
;,  40.1206436,  41.6141205,  43.1032639,  44.3894539,  45.3615913$
;,  45.6144905,  45.0205727,  43.4830284,  41.3165131,  39.0561028$
;,  36.7841721,  34.8829460,  33.2252960,  31.6013317,  30.2279510$
;,  29.0094757,  28.0873394,  27.4200592,  27.0726490,  26.9486694$
;,  27.2879734,  28.2135773,  29.0468178,  26.5301857,  19.4919891$
;,  14.9426184,  12.8888731,  11.8966293,  11.2804928,  10.8809576$
;,  10.5890980,  10.3209915,  10.0815916,   9.8974495,   9.7401304$
;,   9.6096430,   9.5203848,   9.4354925,   9.3487301,   9.2902098$
;,   9.2515221,   9.2242622,   9.2087955,   9.2098846,   9.2018490$
;,   9.1901045,   9.1758842,   9.1766224,   9.2159214,   9.2059631$
;,   9.1819897,   9.1475372,   9.1192532,   9.1246281,   9.1392412$
;,   0.0000000$
; ]
if ( sw_NoPlot eq 0 ) then    oPlot, x_dsp0, tec5, linestyle=4,color=col_min, THICK =line_thick
;
;v6 windX2.5
; tec6=fltarr(nx)
; tec6 =[$
;    0.0000000,   7.5471730,   7.8395762,   8.0856581,   8.4231586$
;,   8.8828840,   9.4369335,  10.0693302,  10.7384844,  11.7319212$
;,  13.0663681,  14.9550982,  17.3666248,  20.3468704,  23.7667770$
;,  27.7144318,  32.4922867,  37.5725670,  43.0288620,  48.6417694$
;,  53.8784714,  58.5305748,  62.7906532,  66.6628113,  69.8208542$
;,  72.2740326,  74.4181290,  75.8253326,  76.2341309,  75.9813461$
;,  74.0552521,  68.0378189,  57.3112450,  45.8223953,  37.6267929$
;,  33.3011665,  31.0489368,  29.9619865,  29.4627666,  29.2875710$
;,  29.3197918,  29.4022217,  29.5507107,  29.8030109,  30.2626514$
;,  30.8401909,  31.5280895,  32.1861382,  32.7330551,  33.1116524$
;,  33.2993507,  33.2783279,  33.1289825,  32.9394684,  32.7738495$
;,  32.6878319,  32.4994926,  32.4607697,  32.5435867,  32.6469345$
;,  32.5480995,  33.0024414,  33.9057045,  35.0026131,  36.3589554$
;,  37.8244247,  39.2167511,  40.6054382,  41.7767334,  42.6738472$
;,  42.8701324,  42.1998444,  40.6537857,  38.4579010,  36.1671982$
;,  33.9294968,  31.9598408,  30.2305317,  28.5303230,  27.0803127$
;,  25.7740154,  24.7096882,  23.8029633,  23.1533356,  22.6667976$
;,  22.5337181,  22.9757538,  24.2421227,  25.6513500,  23.3663177$
;,  17.3309269,  13.5730133,  11.6240635,  10.6832705,  10.1452169$
;,   9.7907238,   9.4997807,   9.2411480,   9.0365620,   8.8696938$
;,   8.7244625,   8.6164484,   8.5206223,   8.4235220,   8.3492203$
;,   8.2946939,   8.2501745,   8.2154560,   8.2002506,   8.1817360$
;,   8.1627102,   8.1457529,   8.1432562,   8.1754675,   8.1772861$
;,   8.1716690,   8.1631508,   8.1581268,   8.1803455,   8.2155561$
;,   0.0000000$
; ]
if ( sw_NoPlot eq 0 ) then oPlot, x_dsp0, tec6, linestyle=2,color=col_min, THICK =line_thick
if ( sw_NoPlot eq 0 ) AND ( sw_output2file eq 1  ) then  output_png, 'Cmptec_'+FileID+'_'+TEST2+'_'+TEST1+'.png'
;
;

if ( sw_NoPlot eq 0 ) then  device, decomposed = 0 ,retain=2
if ( sw_NoPlot eq 0 ) then   window, 3 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

nmf20=nmf2/fac_ne
;print, "nmf20=[$"
;print,(nmf20),format='(5(",",f12.7),"$")'
;print, "]"
if ( sw_NoPlot eq 0 ) then Plot, x_dsp0, nmf20 $
;, yrange=[ 0., 2.5 ], ystyle=1  $ ;20141119
, yrange=[ 0., 2.7 ], ystyle=1  $ ;20141119
;, yrange=[ 0., 10. ], ystyle=1  $
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50., +50. ], xstyle=1  $
, TITLE='nmf2: '+FileID $
,color=text_color, THICK =line_thick, charsize=char_size, charthick=char_thick

;v4: windX1.5
; nmf24=fltarr(nx)
; nmf24 =[$
;   0.0000000,   0.3613207,   0.3689828,   0.3807282,   0.3951636$
;,   0.4149923,   0.4415447,   0.4690003,   0.4975650,   0.5434659$
;,   0.6054812,   0.6781601,   0.7600668,   0.8571669,   0.9529114$
;,   1.0714598,   1.2073280,   1.3252164,   1.4500458,   1.5925921$
;,   1.7083519,   1.8174975,   1.9341639,   2.0361879,   2.1158192$
;,   2.1959784,   2.2744250,   2.3106606,   2.3360472,   2.3782523$
;,   2.3854229,   2.3838415,   2.3996658,   2.3074458,   1.8917722$
;,   1.4999043,   1.3539340,   1.2809116,   1.2284228,   1.1924832$
;,   1.1659143,   1.1478395,   1.1383289,   1.1369349,   1.1471668$
;,   1.1627898,   1.1837015,   1.1936961,   1.1820000,   1.1597043$
;,   1.1190466,   1.0687176,   1.0207462,   0.9830523,   0.9618167$
;,   0.9586117,   0.9630023,   0.9774706,   0.9858128,   0.9895932$
;,   1.0007426,   1.0120174,   1.0609964,   1.1388865,   1.2361945$
;,   1.3446926,   1.4742836,   1.6282512,   1.7806865,   1.9482052$
;,   2.0884964,   2.1934867,   2.2399819,   2.2021401,   2.1527827$
;,   2.0465608,   1.9692335,   1.9161634,   1.8848349,   1.8468266$
;,   1.8024024,   1.7456443,   1.7142593,   1.6928990,   1.6897976$
;,   1.7155622,   1.7828659,   1.7440007,   1.3473557,   1.0297879$
;,   0.8573910,   0.7831838,   0.7604474,   0.7453696,   0.7349336$
;,   0.7246501,   0.7147899,   0.7063718,   0.6974583,   0.6888993$
;,   0.6806874,   0.6723593,   0.6639115,   0.6558794,   0.6488870$
;,   0.6406623,   0.6323183,   0.6288952,   0.6269000,   0.6242282$
;,   0.6226702,   0.6191810,   0.6147550,   0.6131398,   0.6085563$
;,   0.6025261,   0.5930774,   0.5830849,   0.5738415,   0.5691067$
;,   0.0000000$
;]
if ( sw_NoPlot eq 0 ) then oPlot, x_dsp0, nmf24, linestyle=5,color=col_min, THICK =line_thick
;
;v5: windX2
; nmf25=fltarr(nx)
; nmf25 =[$
;   0.0000000,   0.3402883,   0.3466927,   0.3567730,   0.3692365$
;,   0.3855355,   0.4034158,   0.4220643,   0.4428663,   0.4815591$
;,   0.5359309,   0.6044888,   0.6941153,   0.7985827,   0.9187588$
;,   1.0495192,   1.2173563,   1.3665329,   1.5179193,   1.6930141$
;,   1.8323848,   1.9769804,   2.1060200,   2.2200272,   2.3267860$
;,   2.4030938,   2.4746222,   2.5187554,   2.5475581,   2.5639346$
;,   2.5227187,   2.4525611,   2.2526493,   1.8675026,   1.3115402$
;,   1.1443615,   1.0829321,   1.0518079,   1.0281948,   1.0145540$
;,   1.0063421,   1.0025307,   1.0028213,   1.0112643,   1.0300728$
;,   1.0560802,   1.0756820,   1.0814729,   1.0795747,   1.0531021$
;,   1.0080340,   0.9618864,   0.9152631,   0.8831577,   0.8714296$
;,   0.8758361,   0.8881305,   0.9163249,   0.9358673,   0.9383223$
;,   0.9513079,   0.9644731,   1.0090395,   1.0898455,   1.1911013$
;,   1.3133618,   1.4464077,   1.6067660,   1.7862006,   1.9492965$
;,   2.1171763,   2.1768889,   2.2404051,   2.1778548,   2.0494480$
;,   1.9657167,   1.8973218,   1.8556397,   1.7677556,   1.6842963$
;,   1.6128639,   1.5746104,   1.5570916,   1.5373352,   1.5139380$
;,   1.5098345,   1.5398285,   1.6384445,   1.8079718,   1.3167953$
;,   0.9303698,   0.8117903,   0.7632379,   0.7362776,   0.7171435$
;,   0.7016952,   0.6881132,   0.6763145,   0.6651080,   0.6549093$
;,   0.6453529,   0.6371409,   0.6286591,   0.6205831,   0.6143811$
;,   0.6070268,   0.5996935,   0.5936468,   0.5864248,   0.5796317$
;,   0.5729475,   0.5651626,   0.5573758,   0.5496830,   0.5407873$
;,   0.5360639,   0.5316210,   0.5262688,   0.5212577,   0.5197592$
;,   0.0000000$
;]
if ( sw_NoPlot eq 0 ) then oPlot, x_dsp0, nmf25, linestyle=4,color=col_min, THICK =line_thick
;
;v6: windX2.5
; nmf26=fltarr(nx)
; nmf26 =[$
;0.0000000,   0.3338981,   0.3388350,   0.3459444,   0.3542802$
;,   0.3645636,   0.3751982,   0.3863584,   0.4016277,   0.4340949$
;,   0.4741448,   0.5339492,   0.6144755,   0.7218819,   0.8467670$
;,   0.9981155,   1.1794930,   1.3541392,   1.5296842,   1.7311802$
;,   1.8897750,   2.0639365,   2.2003422,   2.3431373,   2.4512258$
;,   2.5203986,   2.6177301,   2.6477928,   2.6503057,   2.6287904$
;,   2.5127559,   2.2993379,   1.7241352,   1.2468541,   1.0092655$
;,   0.9366384,   0.9083075,   0.8970408,   0.8901826,   0.8877230$
;,   0.8897800,   0.8933412,   0.9011157,   0.9130068,   0.9344135$
;,   0.9603309,   0.9786240,   0.9938834,   0.9905694,   0.9627337$
;,   0.9239483,   0.8764210,   0.8312089,   0.8005773,   0.7981445$
;,   0.8097250,   0.8243614,   0.8597949,   0.8869222,   0.8876411$
;,   0.9020225,   0.9166368,   0.9667016,   1.0420190,   1.1501862$
;,   1.2690713,   1.4185511,   1.5771608,   1.7578229,   1.9416622$
;,   2.0757220,   2.1883342,   2.1598520,   2.1368399,   2.0240667$
;,   1.8932600,   1.7874327,   1.7046869,   1.6620679,   1.6116992$
;,   1.5601103,   1.4926331,   1.4468143,   1.3962625,   1.3447546$
;,   1.2979113,   1.2761892,   1.2958323,   1.4412998,   1.6131473$
;,   1.1577018,   0.9151807,   0.7749210,   0.7253411,   0.6919923$
;,   0.6707401,   0.6541993,   0.6397851,   0.6268904,   0.6156965$
;,   0.6053669,   0.5972543,   0.5889362,   0.5810641,   0.5754983$
;,   0.5689368,   0.5625266,   0.5579722,   0.5522734,   0.5469678$
;,   0.5419167,   0.5358627,   0.5299031,   0.5244265,   0.5177225$
;,   0.5108585,   0.5018555,   0.4917645,   0.4816976,   0.4756955$
;,   0.0000000$
;]
if ( sw_NoPlot eq 0 ) then oPlot, x_dsp0, nmf26, linestyle=2,color=col_min, THICK =line_thick
if ( sw_NoPlot eq 0 ) AND ( sw_output2file eq 1  ) then output_png, 'Cmpnmf2_'+FileID+'_'+TEST2+'_'+TEST1+'.png'
;
;---hmf2
if ( sw_NoPlot eq 0 ) then  device, decomposed = 0 ,retain=2
if ( sw_NoPlot eq 0 ) then  window, 4 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

;print, "hmf20=[$"
;print,hmf2,format='(5(",",f12.7),"$")'
;print, "]"
if ( sw_NoPlot eq 0 ) then Plot, x_dsp0, hmf2 $
, yrange=[ 200., 500. ], ystyle=1  $ ;20141119
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50., +50. ], xstyle=1  $
, TITLE='hmf2: '+FileID $
, color=text_color, THICK =line_thick, charsize=char_size, charthick=char_thick

;---
;v4: windX1.5
; hmf24=fltarr(nx)
; hmf24 =[$
;90.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
;, 250.0000000, 250.0000000, 250.0000000, 250.0000000, 250.0000000$
;, 270.0000000, 270.0000000, 290.0000000, 290.0000000, 290.0000000$
;, 310.0000000, 310.0000000, 310.0000000, 330.0000000, 330.0000000$
;, 330.0000000, 350.0000000, 350.0000000, 350.0000000, 350.0000000$
;, 370.0000000, 370.0000000, 370.0000000, 370.0000000, 390.0000000$
;, 390.0000000, 390.0000000, 390.0000000, 410.0000000, 410.0000000$
;, 410.0000000, 390.0000000, 390.0000000, 390.0000000, 390.0000000$
;, 390.0000000, 390.0000000, 390.0000000, 390.0000000, 390.0000000$
;, 410.0000000, 410.0000000, 410.0000000, 410.0000000, 430.0000000$
;, 430.0000000, 430.0000000, 430.0000000, 430.0000000, 410.0000000$
;, 390.0000000, 390.0000000, 350.0000000, 350.0000000, 370.0000000$
;, 370.0000000, 370.0000000, 370.0000000, 370.0000000, 350.0000000$
;, 350.0000000, 330.0000000, 330.0000000, 330.0000000, 310.0000000$
;, 310.0000000, 290.0000000, 290.0000000, 270.0000000, 270.0000000$
;, 270.0000000, 270.0000000, 250.0000000, 250.0000000, 250.0000000$
;, 250.0000000, 250.0000000, 250.0000000, 250.0000000, 250.0000000$
;, 250.0000000, 250.0000000, 230.0000000, 230.0000000, 230.0000000$
;, 230.0000000, 230.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 230.0000000, 230.0000000, 230.0000000$
;, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
;, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
;,  90.0000000$
;]
if ( sw_NoPlot eq 0 ) then oPlot, x_dsp0, hmf24, linestyle=5,color=col_min, THICK =line_thick
;
;---
;v5: windX2.
; hmf25=fltarr(nx)
; hmf25 =[$
; 90.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
;, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 250.0000000$
;, 250.0000000, 270.0000000, 270.0000000, 290.0000000, 290.0000000$
;, 290.0000000, 310.0000000, 310.0000000, 330.0000000, 330.0000000$
;, 330.0000000, 350.0000000, 350.0000000, 370.0000000, 370.0000000$
;, 370.0000000, 390.0000000, 390.0000000, 390.0000000, 390.0000000$
;, 390.0000000, 410.0000000, 410.0000000, 430.0000000, 410.0000000$
;, 410.0000000, 390.0000000, 390.0000000, 390.0000000, 390.0000000$
;, 390.0000000, 390.0000000, 410.0000000, 410.0000000, 410.0000000$
;, 410.0000000, 410.0000000, 410.0000000, 430.0000000, 430.0000000$
;, 430.0000000, 450.0000000, 450.0000000, 430.0000000, 390.0000000$
;, 370.0000000, 370.0000000, 350.0000000, 350.0000000, 370.0000000$
;, 370.0000000, 370.0000000, 370.0000000, 350.0000000, 330.0000000$
;, 330.0000000, 330.0000000, 310.0000000, 310.0000000, 290.0000000$
;, 290.0000000, 290.0000000, 270.0000000, 270.0000000, 250.0000000$
;, 250.0000000, 250.0000000, 250.0000000, 250.0000000, 250.0000000$
;, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
;, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
;,  90.0000000$
;]
if ( sw_NoPlot eq 0 ) then oPlot, x_dsp0, hmf25, linestyle=4,color=col_min, THICK =line_thick
;
;---
;v6: windX2.5
; hmf26=fltarr(nx)
; hmf26 =[$
;  90.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 230.0000000, 230.0000000$
;, 230.0000000, 250.0000000, 250.0000000, 270.0000000, 290.0000000$
;, 290.0000000, 310.0000000, 310.0000000, 330.0000000, 330.0000000$
;, 330.0000000, 350.0000000, 350.0000000, 370.0000000, 370.0000000$
;, 390.0000000, 390.0000000, 390.0000000, 390.0000000, 410.0000000$
;, 410.0000000, 410.0000000, 430.0000000, 430.0000000, 410.0000000$
;, 410.0000000, 410.0000000, 410.0000000, 410.0000000, 410.0000000$
;, 410.0000000, 410.0000000, 410.0000000, 410.0000000, 410.0000000$
;, 410.0000000, 410.0000000, 430.0000000, 430.0000000, 450.0000000$
;, 450.0000000, 450.0000000, 450.0000000, 430.0000000, 390.0000000$
;, 370.0000000, 370.0000000, 350.0000000, 350.0000000, 370.0000000$
;, 370.0000000, 370.0000000, 350.0000000, 330.0000000, 330.0000000$
;, 310.0000000, 310.0000000, 310.0000000, 290.0000000, 290.0000000$
;, 270.0000000, 270.0000000, 270.0000000, 250.0000000, 250.0000000$
;, 250.0000000, 250.0000000, 250.0000000, 230.0000000, 230.0000000$
;, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
;, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
;, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 230.0000000$
;,  90.0000000$
;]
if ( sw_NoPlot eq 0 ) then  oPlot, x_dsp0, hmf26, linestyle=2,color=col_min, THICK =line_thick
if ( sw_NoPlot eq 0 ) AND ( sw_output2file eq 1  ) then output_png, 'Cmphmf2_'+FileID+'_'+TEST2+'_'+TEST1+'.png'
;

if ( sw_debug eq 1 ) then  print,'runID= ',TEST,' ipe_',TEST2,'_',TEST1
if ( sw_debug eq 1 ) then  print, 'time_string=', time_string
if ( sw_SaveTec eq 1 ) then  save, /VARIABLES, filename='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/fig/'+TEST+'/'+TEST1+'_ipe_theia_intel_parallel_'+TEST2+'/TEC_'+TEST+'_ipe_'+TEST2+'_'+TEST1+'_'+time_string+'.sav'

;output to a file
if ( sw_output_nmf2noon eq 1 ) then begin
   if ( n_read eq 0 ) then begin
;   luntmpN=100
      flnmtmpN='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r319/trunk/run/ipe_'+TEST2+'_'+TEST1+'/noon_nmf2.v2.dat'
      openw,luntmpN,flnmtmpN, /GET_LUN
      print, 'noon nmf2 file created:',luntmpN,flnmtmpN
   endif                        ;if n_read eq 0
   printf, luntmpN, mp;mp_output_noonprfl 
print, 'check mp_output=', mp;mp_output_noonprfl 
   printf, luntmpN, uthr, glon_deg[mp]
print, 'ut=', uthr,' glon=', glon_deg[mp]
   printf, luntmpN, lt_hr[mp] ;lt at midpoint
print, 'ltime=', lt_hr[mp] ;lt
   printf, luntmpN, x_dsp0 ;mlat
   printf, luntmpN, nmf2
   printf, luntmpN, hmf2
   printf, luntmpN, tec
endif                  
 

endif ;( sw_plot_tec eq 1 ) then begin

endfor  ;VarType=0,1 do begin  ;0,NPAR-1
endfor ;mp=0,40, 40 do begin



end ; PRO  ctr_lat_ht_tec 
