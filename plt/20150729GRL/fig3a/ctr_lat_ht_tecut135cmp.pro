;20141007 added interpolation for smoother plot and TEC/nmf2/hmf2 calculations 
;TODO
;include parallel V for debug purpose!!!
;20110203: copied from /home/naomi/sandbox/SUPIM/plot_idl/contour_result.pro

 PRO  ctr_lat_ht_tecut135cmp $ ;tmp20160125
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
endif

fac_ne = 1.0E+12
sw_arw_vpara=0
fac_arw_para=0.7

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
  device, decomposed = 0 ,retain=2
  window, 0 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

;t  axis_color = 0.0 ;255.99 $


endif

LOADCT, N_LDCT

;if ( sw_plot_grid eq 1 ) then begin ;20120328
;LOADCT, 0
Plot, xx(0:istop), yy(0:istop)     $
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
Draw_Colorbar, ARY_min0(VarType), ARY_max0(VarType), N_LVLs $
, col_min, col_max, X1, Y1, dX1, dY1, X_SIZE, Y_SIZE, VarType, text_color

;LT
;xyouts, (X0-2.)/X_SIZE, (Y0+dY+2.45)/Y_SIZE, STRTRIM(STRING( zthr, FORMAT='(F5.2)'),1)+'LT'  $
;, charsize=1.0, charthick=1.0, /norm, /noclip


print,'MIN',ary_minz,' MAX',ary_maxz
; add MIN & MAX values
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
help, Tr

d_x=1.0 ;deg                                                                         
d_y=20.;km

trigrid_result=trigrid(x3,y3,z3,tr  $
, [d_x,d_y], [gLATmin,HTmin,gLATmax,HTmax] $
;, NX=12, NY=24  $                                                              
)
size_result=size(trigrid_result)
print, 'check trigrid_result', size_result , min(trigrid_result),max(trigrid_result)

nx=size_result(1) ;111
ny=size_result(2) ;38

x_dsp0=findgen(nx)*d_x + gLATmin
y_dsp0=findgen(ny)*d_y +   HTmin
z_dsp0=(10^(trigrid_result)) /fac_ne
;z_dsp0=trigrid_result ;alog10

zmax=MAX(z_dsp0)
zmin=MIN(z_dsp0)
print, 'Ne^10-12 zmin',zmin,'zmax',zmax
zmax=ARY_max0[VarType]
zmin=ARY_min0[VarType]
  device, decomposed = 0 ,retain=2
  window, 1 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

X_max=+55.
X_min=-X_max
Y_min=90.
Y_max=1000.
char_size = 1.8
char_thick = 2.
line_thick = 3.
text_color=0.

 !P.BACKGROUND=255

contour,z_dsp0,x_dsp0,y_dsp0 $
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
Draw_Colorbar, zmin, zmax, N_LVLs $
, col_min, col_max, X1, Y1, dX1, dY1, X_SIZE, Y_SIZE, VarType, text_color

output_png, 'fig2_'+FileID+'_'+TEST2+'_'+TEST1+'.png'

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

;if n_read eq 48 then begin ;mp=1
if n_read eq 44 then begin  ;mp=41
 print, 'saving at n_read', n_read
 save,filename='tec'+TEST+TEST2+TEST1+'.sav'
endif

;tmp20160116 !P.BACKGROUND=0
  device, decomposed = 0 ,retain=2
  window, 2 ,XSIZE=1000*fac_window,YSIZE=800*fac_window


print,'MAX TEC',MIN(tec/TECU),MAX(tec/TECU)
print,"tec1=[$"
print, tec/TECU, format='(5(",",f12.7),"$")'
print, "]"
Plot, x_dsp0, tec/TECU $
;, yrange=[ 0., 120. ], ystyle=1  $ ;20141119
;, yrange=[ 0., 65. ], ystyle=1  $
, yrange=[ 0., 80. ], ystyle=1  $
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50. , +50. ], xstyle=1  $
, TITLE='TEC: '+FileID $
,color=text_color, THICK =line_thick, charsize=char_size, charthick=char_thick

;v4 windX1.5
 tec4=fltarr(nx)
 tec4 =[$
   0.0000000,  12.7959538,  13.6574163,  14.4672928,  15.5768585$
,  17.1014194,  18.7521076,  20.5170231,  22.0887394,  23.9022369$
,  26.1048584,  28.4244270,  30.9184246,  33.6513557,  36.4258003$
,  39.2388115,  42.3435516,  45.2533035,  48.1234093,  50.9887199$
,  53.6054459,  55.8693886,  57.7820969,  59.6382370,  61.1606636$
,  62.3284340,  63.4630203,  64.2977600,  64.7243271,  65.0658340$
,  65.2824860,  64.9258804,  64.6362381,  64.4207458,  62.4236107$
,  57.2623482,  50.0660820,  44.2846260,  40.3522148,  37.6549187$
,  35.7835159,  34.4695053,  33.6143265,  33.0943146,  32.8746605$
,  32.8886452,  33.1084557,  33.4439049,  33.8468285,  34.2317085$
,  34.5454788,  34.7597275,  34.8811722,  34.9458122,  34.9979439$
,  35.0793419,  34.9569778,  34.9369431,  34.9971962,  35.0236588$
,  34.9409561,  35.4764862,  36.4546356,  37.6295242,  39.0254898$
,  40.4433556,  41.5519218,  42.4413338,  42.9004669,  42.8297005$
,  42.0584068,  40.7883224,  39.2444954,  37.8106308,  36.6794319$
,  35.6559219,  34.8827515,  34.2642708,  33.6992760,  33.3992882$
,  33.3833580,  33.6824875,  34.2621918,  35.0965881,  35.5260925$
,  33.6072121,  26.8184052,  18.8856087,  15.2100563,  13.9605236$
,  12.9973230,  12.4011993,  12.0172949,  11.6833162,  11.4620991$
,  11.2956352,  11.1106529,  10.9477844,  10.8351059,  10.7237234$
,  10.6343651,  10.5835695,  10.5167637,  10.4435253,  10.4026175$
,  10.3627472,  10.3198023,  10.2785625,  10.2429514,  10.1838655$
,  10.1145124,  10.0231819,   9.9387445,   9.8887100,   9.7607803$
,   9.5999880,   9.4163103,   9.2661753,   9.1656199,   9.0911045$
,   0.0000000$
 ]
oPlot, x_dsp0, tec4, linestyle=5 ,color=col_min, THICK =line_thick
;
;v5 windX2
 tec5=fltarr(nx)
 tec5 =[$
   0.0000000,  11.8512850,  12.6410408,  13.4197693,  14.4900961$
,  16.0202312,  17.7482243,  19.6499710,  21.4459076,  23.6126366$
,  26.2754974,  29.1794014,  32.3497963,  35.8707886,  39.4815025$
,  43.1310234,  47.1321983,  50.7804832,  54.2588043,  57.5847321$
,  60.4978294,  62.9253502,  64.7817612,  66.5685959,  67.9509125$
,  68.8938675,  69.7818222,  70.3694458,  70.5281906,  70.4972839$
,  70.3094940,  69.4799728,  68.2575378,  66.5108643,  59.6800270$
,  49.2971878,  41.4702797,  36.8465004,  33.9841156,  32.1241570$
,  31.1091404,  30.5818710,  30.2796803,  30.1387329,  30.2081337$
,  30.4283924,  30.8123531,  31.2783794,  31.7685795,  32.2188683$
,  32.5747910,  32.8053589,  32.9299965,  32.9952049,  33.0514984$
,  33.1446495,  33.0685997,  33.1189461,  33.2834206,  33.4343109$
,  33.4098549,  33.9600906,  34.9265213,  36.0947952,  37.4859238$
,  38.9089394,  40.1022224,  41.1109581,  41.6827354,  41.7197227$
,  40.9306183,  39.4578209,  37.6035728,  35.8156471,  34.4141769$
,  33.2760811,  32.3905754,  31.6663799,  30.9377346,  30.3816185$
,  30.0355740,  30.0910625,  30.4655762,  31.2338734,  32.4800301$
,  34.0860405,  34.1308060,  27.9959621,  17.4372845,  13.4990406$
,  12.0538645,  11.2374916,  10.7445116,  10.3619204,  10.0969915$
,   9.9003878,   9.7074175,   9.5288172,   9.3945951,   9.2752304$
,   9.1705074,   9.0993967,   9.0264864,   8.9467831,   8.8911743$
,   8.8442917,   8.7978773,   8.7536745,   8.7234793,   8.6792612$
,   8.6269798,   8.5640774,   8.5082388,   8.4818468,   8.4037046$
,   8.3019238,   8.1845551,   8.0927420,   8.0392275,   8.0094042$
,   0.0000000$
 ]
oPlot, x_dsp0, tec5, linestyle=4,color=col_min, THICK =line_thick
;
;v6 windX2.5
 tec6=fltarr(nx)
 tec6 =[$
   0.0000000,  10.7760143,  11.4661016,  12.1560717,  13.1091766$
,  14.5118113,  16.1814995,  18.0701866,  19.9830914,  22.3837490$
,  25.3702526,  28.7345581,  32.4679871,  36.6774559,  41.0458298$
,  45.4613228,  50.2732353,  54.5620728,  58.5298271,  62.1731529$
,  65.2431107,  67.7204132,  69.4309616,  71.0849915,  72.2862320$
,  73.0189667,  73.6754761,  74.0374451,  73.9741440,  73.6171570$
,  73.0278625,  71.5784531,  68.5617447,  62.9645538,  51.0070267$
,  40.5389481,  34.4729576,  31.1631279,  29.2499905,  28.2200146$
,  27.8747101,  27.7508163,  27.7251759,  27.7838917,  28.0139751$
,  28.3567104,  28.8402023,  29.3947754,  29.9482155,  30.4483490$
,  30.8381367,  31.0827427,  31.2110634,  31.2774792,  31.3351498$
,  31.4322014,  31.3821564,  31.4639740,  31.6871338,  31.9195156$
,  31.9371910,  32.4899063,  33.4295387,  34.5643616,  35.9261589$
,  37.3275528,  38.5651855,  39.6665535,  40.3596573,  40.5535812$
,  39.8635674,  38.3104324,  36.2799072,  34.2476768,  32.6660004$
,  31.4144611,  30.3911743,  29.5241032,  28.6452312,  27.9063892$
,  27.3150043,  26.9945412,  27.0248508,  27.3994656,  28.1793785$
,  29.6648331,  31.8675537,  32.5233002,  25.8938522,  16.0080166$
,  12.3018894,  10.7974968,  10.0835028,   9.6179523,   9.2995014$
,   9.0670996,   8.8579950,   8.6617775,   8.5075645,   8.3771648$
,   8.2584581,   8.1693125,   8.0865173,   7.9970341,   7.9261971$
,   7.8669271,   7.8102493,   7.7556877,   7.7187004,   7.6740980$
,   7.6236777,   7.5687542,   7.5217957,   7.5011067,   7.4458632$
,   7.3743348,   7.2924213,   7.2303109,   7.1995034,   7.1913705$
,   0.0000000$
 ]
oPlot, x_dsp0, tec6, linestyle=2
;
;

  device, decomposed = 0 ,retain=2
  window, 3 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

print, "nmf21=[$"
print,(nmf2/fac_ne),format='(5(",",f12.7),"$")'
print, "]"
Plot, x_dsp0, nmf2/fac_ne $
;, yrange=[ 0., 2.5 ], ystyle=1  $ ;20141119
, yrange=[ 0., 2.7 ], ystyle=1  $ ;20141119
;, yrange=[ 0., 10. ], ystyle=1  $
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50., +50. ], xstyle=1  $
, TITLE='nmf2: '+FileID

;v4: windX1.5
 nmf24=fltarr(nx)
 nmf24 =[$
   0.0000000,   0.4500801,   0.4734113,   0.5114410,   0.5643194$
,   0.6320712,   0.6949007,   0.7629296,   0.8291726,   0.9079027$
,   0.9906651,   1.0797005,   1.1772623,   1.2816993,   1.3760575$
,   1.4740601,   1.5939497,   1.6850066,   1.7783810,   1.8853980$
,   1.9648321,   2.0534484,   2.1188321,   2.1750588,   2.2274203$
,   2.2699993,   2.3107264,   2.3175511,   2.3305292,   2.3434510$
,   2.3394971,   2.3274910,   2.3303602,   2.3143017,   2.1970782$
,   1.8733760,   1.4986401,   1.3120099,   1.2142328,   1.1491007$
,   1.1014323,   1.0658817,   1.0407476,   1.0241905,   1.0209234$
,   1.0270575,   1.0357206,   1.0434021,   1.0464332,   1.0425822$
,   1.0227679,   0.9942610,   0.9641851,   0.9398462,   0.9227136$
,   0.9197601,   0.9309093,   0.9499072,   0.9634317,   0.9774320$
,   0.9953970,   1.0179372,   1.0674582,   1.1486490,   1.2416546$
,   1.3544842,   1.4773883,   1.6208098,   1.7329565,   1.8444514$
,   1.8753517,   1.8886673,   1.8169787,   1.7904379,   1.7257071$
,   1.6578943,   1.6198456,   1.6299899,   1.6146818,   1.5978957$
,   1.5797049,   1.5634197,   1.5752847,   1.6110976,   1.6751046$
,   1.7422475,   1.6185477,   1.2219845,   0.9472950,   0.8332227$
,   0.7340288,   0.7117263,   0.6991711,   0.6878261,   0.6795449$
,   0.6720169,   0.6643315,   0.6573934,   0.6505306,   0.6434878$
,   0.6365885,   0.6300009,   0.6225404,   0.6152830,   0.6092075$
,   0.6016007,   0.5937954,   0.5861930,   0.5814176,   0.5774364$
,   0.5752937,   0.5709631,   0.5652780,   0.5622215,   0.5562643$
,   0.5488098,   0.5393441,   0.5288777,   0.5184370,   0.5125383$
,   0.0000000$
]
oPlot, x_dsp0, nmf24, linestyle=5
;
;v5: windX2
 nmf25=fltarr(nx)
 nmf25 =[$
   0.0000000,   0.4331189,   0.4559062,   0.4949953,   0.5432344$
,   0.6044158,   0.6719564,   0.7460672,   0.8135659,   0.9017571$
,   1.0077684,   1.1127294,   1.2246655,   1.3594934,   1.4819273$
,   1.6016616,   1.7535434,   1.8602459,   1.9881341,   2.1002934$
,   2.2089400,   2.3052125,   2.3617494,   2.4396160,   2.4770033$
,   2.5073321,   2.5505137,   2.5584736,   2.5528727,   2.5378428$
,   2.5248690,   2.4953854,   2.4528894,   2.3763206,   1.9447047$
,   1.3763175,   1.1515306,   1.0565464,   1.0019145,   0.9662305$
,   0.9429864,   0.9265958,   0.9149456,   0.9093772,   0.9124155$
,   0.9238915,   0.9399556,   0.9546750,   0.9608330,   0.9508922$
,   0.9313244,   0.9039539,   0.8723505,   0.8446351,   0.8291925$
,   0.8335558,   0.8492194,   0.8813850,   0.9101983,   0.9253331$
,   0.9469061,   0.9689821,   1.0195030,   1.0971211,   1.1972743$
,   1.3211892,   1.4569975,   1.6049148,   1.7569777,   1.8729113$
,   1.9269444,   1.9076660,   1.8454360,   1.7151328,   1.6507666$
,   1.5944048,   1.5658144,   1.5517155,   1.5071017,   1.4628161$
,   1.4491259,   1.4321569,   1.4354534,   1.4495045,   1.4824188$
,   1.5626097,   1.7074443,   1.7851341,   1.1911956,   0.8249844$
,   0.7331010,   0.6880270,   0.6667898,   0.6505077,   0.6383312$
,   0.6286798,   0.6195782,   0.6112953,   0.6037208,   0.5963330$
,   0.5892531,   0.5834807,   0.5767635,   0.5701962,   0.5653463$
,   0.5590392,   0.5525646,   0.5471324,   0.5405784,   0.5342214$
,   0.5283183,   0.5213321,   0.5141131,   0.5069962,   0.4986993$
,   0.4902949,   0.4802265,   0.4702667,   0.4643295,   0.4618016$
,   0.0000000$
]
oPlot, x_dsp0, nmf25, linestyle=4
;
;v6: windX2.5
 nmf26=fltarr(nx)
 nmf26 =[$
   0.0000000,   0.4112466,   0.4291663,   0.4593571,   0.4967637$
,   0.5575349,   0.6207044,   0.6845172,   0.7599383,   0.8565385$
,   0.9642867,   1.0917258,   1.2202389,   1.3732895,   1.5221817$
,   1.6649719,   1.8449740,   1.9781443,   2.1190183,   2.2538514$
,   2.3664060,   2.4682183,   2.5349967,   2.5977852,   2.6455493$
,   2.6721323,   2.6890793,   2.6985435,   2.6990771,   2.6740282$
,   2.6286981,   2.5620885,   2.4084768,   2.2113230,   1.4096917$
,   1.0419916,   0.9198166,   0.8720574,   0.8474348,   0.8331075$
,   0.8258255,   0.8215153,   0.8191248,   0.8202768,   0.8288721$
,   0.8431489,   0.8603520,   0.8748501,   0.8826458,   0.8827583$
,   0.8616282,   0.8301873,   0.7970970,   0.7688310,   0.7525352$
,   0.7625778,   0.7846332,   0.8238764,   0.8598704,   0.8731382$
,   0.8968320,   0.9211689,   0.9732555,   1.0518887,   1.1577297$
,   1.2786829,   1.4245290,   1.5814855,   1.7406299,   1.8791484$
,   1.9478117,   1.9322877,   1.8105066,   1.7338510,   1.6074367$
,   1.5237033,   1.4754754,   1.4661252,   1.4474498,   1.4197353$
,   1.3891947,   1.3522137,   1.3301231,   1.3101643,   1.3021420$
,   1.3225070,   1.4059974,   1.6082876,   1.7398655,   1.0153364$
,   0.7800018,   0.6715744,   0.6348338,   0.6132408,   0.5970812$
,   0.5853004,   0.5751623,   0.5659040,   0.5576843,   0.5501068$
,   0.5429649,   0.5375693,   0.5313911,   0.5253762,   0.5212462$
,   0.5158580,   0.5103685,   0.5063465,   0.5012894,   0.4963831$
,   0.4919255,   0.4864962,   0.4808891,   0.4758596,   0.4696660$
,   0.4631726,   0.4553839,   0.4462891,   0.4366849,   0.4298670$
,   0.0000000$
]
oPlot, x_dsp0, nmf26, linestyle=2
;
;---hmf2
  device, decomposed = 0 ,retain=2
  window, 4 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

print, "hmf21=[$"
print,hmf2,format='(5(",",f12.7),"$")'
print, "]"
Plot, x_dsp0, hmf2 $
, yrange=[ 200., 500. ], ystyle=1  $ ;20141119
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50., +50. ], xstyle=1  $
, TITLE='hmf2: '+FileID
;---
;v4: windX1.5
 hmf24=fltarr(nx)
 hmf24 =[$
  90.0000000, 250.0000000, 250.0000000, 250.0000000, 270.0000000$
, 270.0000000, 270.0000000, 290.0000000, 290.0000000, 290.0000000$
, 290.0000000, 310.0000000, 310.0000000, 310.0000000, 310.0000000$
, 330.0000000, 330.0000000, 330.0000000, 350.0000000, 350.0000000$
, 350.0000000, 370.0000000, 370.0000000, 370.0000000, 390.0000000$
, 390.0000000, 390.0000000, 390.0000000, 410.0000000, 410.0000000$
, 410.0000000, 410.0000000, 410.0000000, 430.0000000, 430.0000000$
, 430.0000000, 430.0000000, 410.0000000, 410.0000000, 410.0000000$
, 410.0000000, 410.0000000, 410.0000000, 410.0000000, 410.0000000$
, 410.0000000, 410.0000000, 410.0000000, 430.0000000, 430.0000000$
, 430.0000000, 430.0000000, 430.0000000, 430.0000000, 410.0000000$
, 410.0000000, 390.0000000, 390.0000000, 390.0000000, 370.0000000$
, 370.0000000, 390.0000000, 370.0000000, 370.0000000, 350.0000000$
, 350.0000000, 330.0000000, 330.0000000, 310.0000000, 310.0000000$
, 290.0000000, 290.0000000, 290.0000000, 270.0000000, 270.0000000$
, 270.0000000, 270.0000000, 250.0000000, 250.0000000, 250.0000000$
, 250.0000000, 250.0000000, 250.0000000, 250.0000000, 250.0000000$
, 250.0000000, 250.0000000, 230.0000000, 230.0000000, 230.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 230.0000000, 230.0000000$
, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
,  90.0000000$
]
oPlot, x_dsp0, hmf24, linestyle=5
;
;---
;v5: windX2.
 hmf25=fltarr(nx)
 hmf25 =[$
  90.0000000, 250.0000000, 250.0000000, 250.0000000, 250.0000000$
, 250.0000000, 270.0000000, 270.0000000, 270.0000000, 290.0000000$
, 290.0000000, 290.0000000, 310.0000000, 310.0000000, 310.0000000$
, 330.0000000, 330.0000000, 330.0000000, 350.0000000, 350.0000000$
, 370.0000000, 370.0000000, 390.0000000, 390.0000000, 390.0000000$
, 410.0000000, 410.0000000, 410.0000000, 410.0000000, 430.0000000$
, 430.0000000, 430.0000000, 430.0000000, 430.0000000, 450.0000000$
, 450.0000000, 430.0000000, 430.0000000, 430.0000000, 430.0000000$
, 430.0000000, 430.0000000, 430.0000000, 430.0000000, 430.0000000$
, 430.0000000, 430.0000000, 430.0000000, 430.0000000, 430.0000000$
, 450.0000000, 450.0000000, 450.0000000, 430.0000000, 410.0000000$
, 390.0000000, 390.0000000, 350.0000000, 350.0000000, 370.0000000$
, 370.0000000, 370.0000000, 370.0000000, 350.0000000, 350.0000000$
, 330.0000000, 330.0000000, 310.0000000, 310.0000000, 290.0000000$
, 290.0000000, 270.0000000, 270.0000000, 270.0000000, 250.0000000$
, 250.0000000, 250.0000000, 250.0000000, 250.0000000, 250.0000000$
, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 230.0000000, 230.0000000, 230.0000000$
,  90.0000000$
]
oPlot, x_dsp0, hmf25, linestyle=4
;
;---
;v6: windX2.5
 hmf26=fltarr(nx)
 hmf26 =[$
  90.0000000, 230.0000000, 230.0000000, 230.0000000, 250.0000000$
, 250.0000000, 250.0000000, 250.0000000, 270.0000000, 270.0000000$
, 290.0000000, 290.0000000, 290.0000000, 310.0000000, 310.0000000$
, 330.0000000, 330.0000000, 350.0000000, 350.0000000, 370.0000000$
, 370.0000000, 390.0000000, 390.0000000, 410.0000000, 410.0000000$
, 410.0000000, 410.0000000, 430.0000000, 430.0000000, 430.0000000$
, 430.0000000, 430.0000000, 450.0000000, 450.0000000, 470.0000000$
, 450.0000000, 430.0000000, 430.0000000, 430.0000000, 430.0000000$
, 430.0000000, 430.0000000, 430.0000000, 430.0000000, 430.0000000$
, 430.0000000, 430.0000000, 430.0000000, 450.0000000, 450.0000000$
, 450.0000000, 450.0000000, 450.0000000, 450.0000000, 430.0000000$
, 390.0000000, 370.0000000, 350.0000000, 350.0000000, 370.0000000$
, 370.0000000, 370.0000000, 350.0000000, 350.0000000, 330.0000000$
, 330.0000000, 310.0000000, 310.0000000, 290.0000000, 290.0000000$
, 270.0000000, 270.0000000, 250.0000000, 250.0000000, 250.0000000$
, 250.0000000, 250.0000000, 230.0000000, 230.0000000, 230.0000000$
, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 230.0000000$
, 230.0000000, 230.0000000, 230.0000000, 230.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
, 210.0000000, 210.0000000, 210.0000000, 210.0000000, 210.0000000$
,  90.0000000$
]
oPlot, x_dsp0, hmf26, linestyle=2
;

print,'runID= ',TEST,' ipe_',TEST2,'_',TEST1
print, 'time_string=', time_string
save, /VARIABLES, filename=TEST+'_ipe_'+TEST2+'_'+TEST1+'_'+time_string+'.sav'

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
