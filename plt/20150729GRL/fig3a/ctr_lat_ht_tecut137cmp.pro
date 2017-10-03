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
print,"tec1=[$"
print, tec/TECU, format='(5(",",f12.7),"$")'
print, "]"
Plot, x_dsp0, tec/TECU $
;, yrange=[ 0., 120. ], ystyle=1  $ ;20141119
, yrange=[ 0., 65. ], ystyle=1  $
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50. , +50. ], xstyle=1  $
, TITLE='TEC: '+FileID

;v4 windX1.5
 tec4=fltarr(nx)
 tec4 =[$
   0.0000000,   8.4231319,   8.7300930,   8.9576454,   9.2582464$
,   9.6260843,  10.0422678,  10.4876366,  10.9023256,  11.4622946$
,  12.1896057,  13.1263371,  14.3019390,  15.6730137,  17.3770218$
,  19.4336205,  22.0562553,  24.9485950,  28.1930294,  31.6946831$
,  35.1857452,  38.6553879,  42.3835411,  46.0913620,  49.6024170$
,  52.9693718,  56.0814934,  58.5961342,  60.6776886,  62.5895271$
,  63.6935654,  63.6781120,  62.8528862,  57.7532959,  50.3017960$
,  44.9486732,  41.8866882,  40.0244751,  38.8215714,  37.9993706$
,  37.5961761,  37.3902473,  37.4468613,  37.7040329,  38.2428017$
,  38.9693451,  39.7029762,  40.3481369,  40.7992325,  40.9849091$
,  40.9169693,  40.5968285,  40.1485062,  39.6794548,  39.2646637$
,  38.9798164,  38.6167030,  38.4156227,  38.2831230,  38.1715546$
,  37.9933701,  38.4525185,  39.4937744,  40.7168007,  42.2550163$
,  43.9952087,  45.6845245,  47.4595451,  49.1582222,  50.6554031$
,  51.6569061,  51.9522247,  51.2084427,  49.4933891,  47.0751038$
,  44.0419159,  41.3329391,  38.8222389,  36.4489136,  34.3796959$
,  32.5463638,  31.0215454,  29.7076969,  28.6691341,  27.7652283$
,  27.1631680,  26.7967033,  25.9484673,  22.7060242,  19.1223316$
,  16.5787163,  15.1704235,  14.3044243,  13.6472921,  13.2152386$
,  12.8718472,  12.5336723,  12.2421494,  12.0158386,  11.8181515$
,  11.6741476,  11.5769920,  11.4778204,  11.3823605,  11.3216944$
,  11.2922392,  11.2875910,  11.2993450,  11.3233242,  11.3398886$
,  11.3561392,  11.3748598,  11.4138336,  11.4988976,  11.5235920$
,  11.5256414,  11.5179377,  11.4918232,  11.4959745,  11.4892664$
,   0.0000000$
 ]
oPlot, x_dsp0, tec4, linestyle=5
;
;v5 windX2
 tec5=fltarr(nx)
 tec5 =[$
   0.0000000,   7.3158884,   7.5404043,   7.7186489,   7.9528327$
,   8.2414484,   8.5751257,   8.9396257,   9.2980480,   9.7814608$
,  10.4075928,  11.2133579,  12.2302160,  13.4378719,  15.0221233$
,  17.1087208,  19.9061584,  23.1484070,  26.8425732,  30.8697701$
,  34.9406662,  39.0660210,  43.4825172,  47.8460312,  51.9634361$
,  55.8534470,  59.3421249,  62.0184479,  63.8971252,  65.0013733$
,  63.9777412,  59.3029137,  52.1554565,  44.3061905,  39.0264053$
,  36.2933884,  34.8371849,  34.0400810,  33.6350441,  33.4606552$
,  33.5090065,  33.6220856,  33.9016380,  34.3280869,  34.9991760$
,  35.7964516,  36.5790062,  37.2275352,  37.6617432,  37.8367500$
,  37.7784538,  37.4899101,  37.0914879,  36.6852608,  36.3382530$
,  36.1131516,  35.8186646,  35.6899300,  35.6514511,  35.6484070$
,  35.4727097,  35.8854675,  36.8294296,  37.9578438,  39.3729591$
,  40.9482574,  42.4735069,  44.0516586,  45.4980011,  46.7460060$
,  47.4453316,  47.4267807,  46.4414749,  44.5637665,  42.1392136$
,  39.2338638,  36.5643349,  34.1113052,  31.7228756,  29.6581402$
,  27.8010788,  26.2231884,  24.7862396,  23.5832729,  22.4761181$
,  21.6056194,  21.0947437,  21.1733570,  21.7725315,  20.3901958$
,  17.0499954,  14.7424545,  13.3222761,  12.4465008,  11.9074020$
,  11.5158968,  11.1686277,  10.8650055,  10.6227570,  10.4205914$
,  10.2581215,  10.1387472,  10.0293818,   9.9230270,   9.8444643$
,   9.7946911,   9.7640638,   9.7484512,   9.7535181,   9.7595491$
,   9.7733011,   9.7956762,   9.8374662,   9.9184446,   9.9671459$
,  10.0040951,  10.0376930,  10.0562716,  10.0965481,  10.1328163$
,   0.0000000$
 ]
oPlot, x_dsp0, tec5, linestyle=4
;
;v6 windX2.5
 tec6=fltarr(nx)
 tec6 =[$
   0.0000000,   6.5372210,   6.7053804,   6.8474054,   7.0327044$
,   7.2607899,   7.5267391,   7.8233275,   8.1306086,   8.5439405$
,   9.0789089,   9.7702837,  10.6349688,  11.6735964,  13.0369129$
,  14.9704447,  17.6082649,  20.9475231,  24.8545895,  29.1677227$
,  33.5879784,  38.1609573,  43.0673523,  47.9211998,  52.5049553$
,  56.7768898,  60.4657707,  63.0153122,  64.0950775,  63.3968163$
,  59.2255211,  50.9789925,  41.9747772,  35.7633133,  32.6765022$
,  31.2639980,  30.5760365,  30.3036728,  30.2539062,  30.3317280$
,  30.5520020,  30.7979755,  31.1670685,  31.6608868,  32.3860550$
,  33.2003555,  33.9954872,  34.6355705,  35.0542068,  35.2241745$
,  35.1768761,  34.9170876,  34.5594788,  34.2006454,  33.9009781$
,  33.7122841,  33.4568100,  33.3617516,  33.3817902,  33.4556313$
,  33.2866020,  33.6597404,  34.5088730,  35.5381126,  36.8256340$
,  38.2296143,  39.5706558,  40.9311600,  42.0945702,  43.0431633$
,  43.4179420,  43.0734215,  41.9400291,  40.0889435,  37.8564644$
,  35.2618179,  32.7477531,  30.3980083,  28.0629959,  26.0468712$
,  24.2251835,  22.6746330,  21.2344265,  19.9708252,  18.7680283$
,  17.7118912,  16.8858051,  16.3679142,  16.3650475,  16.4920292$
,  15.2088604,  13.7965612,  12.5971012,  11.7149601,  11.1222591$
,  10.6909523,  10.3320341,  10.0170689,   9.7613134,   9.5545444$
,   9.3806868,   9.2469616,   9.1302233,   9.0164127,   8.9251995$
,   8.8590183,   8.8073511,   8.7676430,   8.7500963,   8.7367830$
,   8.7339792,   8.7415342,   8.7679548,   8.8301983,   8.8777008$
,   8.9213982,   8.9675875,   9.0021362,   9.0545082,   9.1127892$
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
, yrange=[ 0., 2.5 ], ystyle=1  $ ;20141119
;, yrange=[ 0., 10. ], ystyle=1  $
;, xrange=[ gLATMIN, gLATMAX ], xstyle=1  $
, xrange=[ -50., +50. ], xstyle=1  $
, TITLE='nmf2: '+FileID

;v4: windX1.5
 nmf24=fltarr(nx)
 nmf24 =[$
   0.0000000,   0.3412536,   0.3464689,   0.3537920,   0.3620267$
,   0.3715382,   0.3818711,   0.3916064,   0.4011605,   0.4166085$
,   0.4347579,   0.4636621,   0.5056831,   0.5530950,   0.6194177$
,   0.7040072,   0.8117992,   0.9294876,   1.0516531,   1.1932480$
,   1.3189977,   1.4458995,   1.5804197,   1.6984619,   1.8239344$
,   1.9285319,   2.0226493,   2.0901988,   2.1615982,   2.2334199$
,   2.2364531,   2.2308857,   2.2210512,   1.9037406,   1.5451666$
,   1.4107879,   1.3542730,   1.3199248,   1.2939528,   1.2784107$
,   1.2704097,   1.2693167,   1.2777456,   1.2923769,   1.3146805$
,   1.3421762,   1.3596903,   1.3482140,   1.3272493,   1.2698269$
,   1.2068323,   1.1395718,   1.0816011,   1.0442640,   1.0300955$
,   1.0257548,   1.0285020,   1.0445621,   1.0503885,   1.0456705$
,   1.0529981,   1.0603769,   1.1054450,   1.1871036,   1.2828482$
,   1.3903911,   1.5205815,   1.6671672,   1.8163085,   1.9948514$
,   2.1740704,   2.3119800,   2.4610040,   2.4760518,   2.4402163$
,   2.3449798,   2.2494016,   2.1686301,   2.0379415,   1.9624259$
,   1.9146577,   1.8409028,   1.7924718,   1.7405918,   1.6961101$
,   1.6591398,   1.6457140,   1.6729689,   1.5640962,   1.2542053$
,   1.0182054,   0.9166630,   0.8662610,   0.8299577,   0.8054618$
,   0.7865152,   0.7680177,   0.7530956,   0.7405974,   0.7285945$
,   0.7185590,   0.7124316,   0.7055300,   0.6990438,   0.6956626$
,   0.6914880,   0.6871875,   0.6866746,   0.6840844,   0.6812048$
,   0.6787292,   0.6745421,   0.6697805,   0.6676894,   0.6619898$
,   0.6541533,   0.6461314,   0.6362195,   0.6264459,   0.6206281$
,   0.0000000$
]
oPlot, x_dsp0, nmf24, linestyle=5
;
;v5: windX2
 nmf25=fltarr(nx)
 nmf25 =[$
   0.0000000,   0.3320023,   0.3360732,   0.3413420,   0.3470027$
,   0.3532463,   0.3598708,   0.3661111,   0.3724335,   0.3816238$
,   0.3959732,   0.4166423,   0.4449670,   0.4793352,   0.5363194$
,   0.6144809,   0.7251627,   0.8466945,   0.9919298,   1.1502012$
,   1.2867134,   1.4390134,   1.5992857,   1.7351764,   1.8804487$
,   2.0000710,   2.1057796,   2.1764827,   2.2324736,   2.2771385$
,   2.1653752,   1.9535129,   1.6450890,   1.2983902,   1.1865277$
,   1.1461722,   1.1261140,   1.1148916,   1.1078435,   1.1095941$
,   1.1160550,   1.1253451,   1.1411821,   1.1607566,   1.1860223$
,   1.2175944,   1.2321591,   1.2250904,   1.1969010,   1.1434948$
,   1.0833317,   1.0214257,   0.9738178,   0.9503995,   0.9451590$
,   0.9476820,   0.9542663,   0.9812988,   0.9977670,   0.9913397$
,   0.9994775,   1.0076799,   1.0551715,   1.1330916,   1.2362759$
,   1.3529593,   1.4777160,   1.6393768,   1.8081068,   1.9693416$
,   2.1541777,   2.2761238,   2.3672531,   2.4113770,   2.3140781$
,   2.1575813,   2.0926466,   2.0465744,   1.9443716,   1.8456669$
,   1.7449698,   1.6443017,   1.5594046,   1.5258361,   1.4884565$
,   1.4539932,   1.4244726,   1.4110049,   1.4793978,   1.4803467$
,   1.1293159,   0.9390953,   0.8520587,   0.8157932,   0.7927592$
,   0.7739874,   0.7578619,   0.7445272,   0.7302684,   0.7179439$
,   0.7064369,   0.6950152,   0.6845503,   0.6746280,   0.6658635$
,   0.6565638,   0.6473627,   0.6390704,   0.6297644,   0.6210466$
,   0.6120475,   0.6023172,   0.6007835,   0.6016258,   0.5996724$
,   0.5960921,   0.5926777,   0.5870226,   0.5811087,   0.5784110$
,   0.0000000$
]
oPlot, x_dsp0, nmf25, linestyle=4
;
;v6: windX2.5
 nmf26=fltarr(nx)
 nmf26 =[$
   0.0000000,   0.3227344,   0.3262072,   0.3308451,   0.3358091$
,   0.3412461,   0.3471886,   0.3528912,   0.3588189,   0.3674965$
,   0.3774417,   0.3897918,   0.4062223,   0.4313469,   0.4733759$
,   0.5378550,   0.6337946,   0.7590593,   0.8985214,   1.0718293$
,   1.2193897,   1.3804369,   1.5590109,   1.7135366,   1.8694022$
,   1.9991138,   2.1101122,   2.1718659,   2.1845891,   2.1665308$
,   1.9012173,   1.5707839,   1.1922293,   1.0217779,   0.9838189$
,   0.9712743,   0.9669224,   0.9682929,   0.9715281,   0.9806547$
,   0.9937350,   1.0080340,   1.0274523,   1.0494685,   1.0761901$
,   1.1105691,   1.1227608,   1.1192652,   1.0869131,   1.0397471$
,   0.9809412,   0.9243055,   0.8856571,   0.8735159,   0.8762305$
,   0.8808413,   0.8908416,   0.9229770,   0.9466067,   0.9384645$
,   0.9473783,   0.9563767,   1.0092566,   1.0837981,   1.1872653$
,   1.3049902,   1.4467286,   1.5942507,   1.7685089,   1.9413297$
,   2.0653419,   2.2068503,   2.2492385,   2.2145593,   2.1899252$
,   2.0978780,   1.9978632,   1.9140776,   1.7751426,   1.6482891$
,   1.6090434,   1.5434608,   1.4982747,   1.4414581,   1.3717372$
,   1.2968458,   1.2275652,   1.1643612,   1.1441730,   1.1844581$
,   1.0650979,   0.9886418,   0.8957259,   0.8376918,   0.7989599$
,   0.7706277,   0.7487161,   0.7304047,   0.7125442,   0.6977417$
,   0.6841580,   0.6720169,   0.6608845,   0.6504263,   0.6418663$
,   0.6328795,   0.6242172,   0.6171395,   0.6090136,   0.6015267$
,   0.5940693,   0.5855125,   0.5771664,   0.5692880,   0.5597565$
,   0.5495792,   0.5439637,   0.5413469,   0.5383111,   0.5378916$
,   0.0000000$
]
oPlot, x_dsp0, nmf26, linestyle=2
;

;;;hmf2 is not available!!!
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
