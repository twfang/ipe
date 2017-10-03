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
, COLOR=text_color ;$

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

 !P.BACKGROUND=0
  device, decomposed = 0 ,retain=2
  window, 2 ,XSIZE=1000*fac_window,YSIZE=800*fac_window


print,'MAX TEC',MIN(tec/TECU),MAX(tec/TECU)
print,' tec/TECU=', tec/TECU
Plot, x_dsp0, tec/TECU $
;, yrange=[ 0., 120. ], ystyle=1  $ ;20141119
, yrange=[ 0., 65. ], ystyle=1  $
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50. , +50. ], xstyle=1  $
, TITLE='TEC: '+FileID

;add 80_378(windX1.5 jid79)
 tec1=fltarr(nx)
 tec1 =[  9.80000e-11,      8.39871,      8.69310 ,     8.91015 $
     , 9.19831,      9.55560,      9.95644,      10.3830,      10.7837 $
     , 11.3158,      12.0010,      12.8770,      13.9709,      15.2434$
     , 16.8373,      18.7900,      21.3054,      24.0953,      27.2247$
     , 30.6070,      33.9891,      37.3767,      41.0385,      44.6938$
     , 48.1804,      51.5467,      54.6643,      57.2036,      59.3464$
     , 61.3238,      62.5078,      62.6668,      62.0239,      57.1061$
     , 49.8870,      44.6828,      41.6975,      39.8858,      38.7253$
     , 37.9316,      37.5493,      37.3625,      37.4405,      37.7158$
     , 38.2725,      39.0121,      39.7503,      40.3974,      40.8469$
     , 41.0298,      40.9599,      40.6383,      40.1890,      39.7193$
     , 39.3036,      39.0176,      38.6529,      38.4494,      38.3136$
     , 38.1994,      38.0183,      38.4731,      39.5118,      40.7305$
     , 42.2646,      44.0007,      45.6833,      47.4506,      49.1393$
     , 50.6244,      51.6151,      51.8970,      51.1391,      49.4124$
     , 46.9823,      43.9389,      41.2277,      38.7142,      36.3386$
     , 34.2583,      32.4083,      30.8587,      29.5137,      28.4387$
     , 27.4956,      26.8516,      26.4403,      25.5758,      22.4648$
     , 19.0862,      16.6328,      15.2602,      14.4041,      13.7458$
     , 13.3113,      12.9640,      12.6206,      12.3252,      12.0955$
     , 11.8941,      11.7475,      11.6480,      11.5463,      11.4486$
     , 11.3862,      11.3554,      11.3496,      11.3605,      11.3846$
     , 11.4021,      11.4213,      11.4440,      11.4865,      11.5727$
     , 11.6006,      11.6063,      11.6007,      11.5802,      11.5911$
     , 11.5894,  9.80000e-11 ]
oPlot, x_dsp0, tec1, linestyle=5

  device, decomposed = 0 ,retain=2
  window, 3 ,XSIZE=1000*fac_window,YSIZE=800*fac_window

print, 'nmf2/fac_ne=', nmf2/fac_ne
Plot, x_dsp0, nmf2/fac_ne $
, yrange=[ 0., 2.5 ], ystyle=1  $ ;20141119
;, yrange=[ 0., 10. ], ystyle=1  $
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50., +50. ], xstyle=1  $
, TITLE='nmf2: '+FileID

;add 80_378(windX1.5 jid79)
 nmf21=fltarr(nx)
 nmf21 =[ 1.00000e-12,     0.341161,     0.346208 ,    0.353282 $
   ,  0.361261  ,   0.370534,     0.380338,     0.389636,     0.398821$
   ,  0.413207  ,   0.430026,     0.454500,     0.493296,     0.537677$
   ,  0.596334  ,   0.678687,     0.778834,     0.894199,      1.01443$
   ,   1.14900 ,     1.27330,      1.39962,      1.53462,      1.64892$
   ,   1.77276 ,     1.88024,      1.97771,      2.04543,      2.11369$
   ,   2.19021 ,     2.19909,      2.20309,      2.20533,      1.89030$
   ,   1.53853 ,     1.41300,      1.35754,      1.32087,      1.29282$
   ,   1.27865 ,     1.27213,      1.27241,      1.28258,      1.29860$
   ,   1.32245 ,     1.34671,      1.36486,      1.35085,      1.33010$
   ,   1.27282 ,     1.20818,      1.14135,      1.08492,      1.04819$
   ,   1.03372 ,     1.02907,      1.03132,      1.04670,      1.05202$
   ,   1.04704 ,     1.05400,      1.06102,      1.10594,      1.18712$
   ,   1.28257 ,     1.38847,      1.51804,      1.66440,      1.81379$
   ,   1.98816 ,     2.16907,      2.30018,      2.45239,      2.47274$
   ,   2.42617 ,     2.33437,      2.24230,      2.16567,      2.03890$
    ,  1.94966 ,     1.90620,      1.83478,      1.78831,      1.73709$
   ,   1.69244 ,     1.65465,      1.63907,      1.66152,      1.53931$
   ,   1.24793 ,     1.02036,     0.922626,     0.873478,     0.837307$
   ,  0.812841 ,    0.793731,     0.774970,     0.759858,     0.747144$
   ,  0.734901,     0.724655,     0.718299,     0.711205,     0.704544$
    , 0.701008,     0.696685,     0.692239,     0.691579,     0.688846$
    , 0.685841,     0.683328,     0.679154,     0.674478,     0.672253$
    , 0.666557,     0.658878,     0.650776,     0.640844,     0.631077$
    , 0.625231,  1.00000e-12]
oPlot, x_dsp0, nmf21, linestyle=5
print,'runID ',TEST,' ipe_',TEST2,'_',TEST1

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
