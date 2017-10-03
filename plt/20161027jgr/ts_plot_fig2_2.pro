;070306: modified to include P_dyn --> Fig2_2
;052506: modified to include Fig3 --> Fig2_1
;052306: adela drift with the neural network is replaced with the slope version
pro ts_plot_fig2_2
 compile_opt strictarr

DIR_DISP = '/Users/naomi/jastp_special2006/figure/finals062906/'


;*** determine the thickness of arrow line
 !P.thick=4.5
; choose fonts: 1:TrueType Fonts, -1:Vector-drawn, 0:graphic device specified
!P.FONT=1
char_size = 1.5
char_thick = 2.0
;choose color
 col_max = 256.0
 col_min =   0.0
N_LDCT = 39L   ;rainbow+black
loadct, N_LDCT 


;*** parameters for ps plot
FILE_DISP_plot = DIR_DISP+'test_fig2s_2.ps'
; window size
 X_SIZE=29.7
 Y_SIZE=20.5 ;1.0
 X0=8.0
 dX=15.0
;****

SET_PLOT, 'ps'
;;;;;;;;;;;; IF(DEVICE_TYPE eq 1) then $
 DEVICE, BITS_PER_PIXEL=24,/COLOR, $
 FILENAME=FILE_DISP_plot,       /LANDSCAPE, $
 SET_FONT='Courier',       /TT_FONT,   $
 XSIZE=X_SIZE, YSIZE=Y_SIZE, XOFFSET=0.0, YOFFSET=X_SIZE


;DEVICE, RETAIN=2, DECOMPOSED=0
;WINDOW,0,XSIZE=900,YSIZE=700
;!p.multi=[0,1,3]



;0) readin rcm PCP_kV  ...
restore, filename='rcmpcp_inp033101.sav'
help, /structure, RCMPCP_INP_ASCII
size_result = SIZE( RCMPCP_INP_ASCII.UT_SEC )
jmax = size_result[1]  ;=718
print, 'jmax', jmax
;dbg: print, min( RCMPCP_INP_ASCII.UT_SEC )/3600., max( RCMPCP_INP_ASCII.UT_SEC )/3600.
jmin = 347L
print, 'min',min(RCMPCP_INP_ASCII.PCP_KV[jmin:jmax-1]), ' max',max(RCMPCP_INP_ASCII.PCP_KV[jmin:jmax-1])
print, jmin,RCMPCP_INP_ASCII.UT_SEC[jmin-1:jmin+1]/3600.


;1st panel: (1) left: PCP[kV]
Y00 = 2.50
dY = 5.30;14.0
dY_margin = 0.99
Y0 = Y00 +(dY +dY_margin)*2.0


Y_MIN=  +0.00 
Y_MAX=+250.00 
X_MIN=  0.
X_MAX=+24.
Main_Title = 'Figure 2s_2:   Jicamarca LON=283[deg.]'

plot, (RCMPCP_INP_ASCII.UT_SEC[jmin:jmax-1]/3600. -24.), RCMPCP_INP_ASCII.PCP_KV[jmin:jmax-1] $
, yrange=[ Y_MIN, Y_MAX ], ystyle=8  $
, xrange=[ X_MIN, X_MAX ], xstyle=1  $
;, XTITLE='UT [hours]'  $
, YTITLE='PCP [kV]'  $
, charsize=char_size, charthick=char_thick  $
, LINESTYLE=0    $
;, XTICKINTERVAL=3.0  $
, pos=[X0/X_SIZE, Y0/Y_SIZE, (X0+dX)/X_SIZE, (Y0+dY)/Y_SIZE]  $
, title=Main_Title  $
, XTICKs=8, xminor=3  $
, XTICKLEN=1, XGRIDSTYLE=1  $ 
, YTICKs=5, Yminor=5        $ ;, xtickname=['0',' ','50',' ','100',' ','150',' ','200',' ','250' ]  $
, YTICKLEN=1, YGRIDSTYLE=1 ;$


;1) readin rcm inputs, IMF Bz ...
restore, filename='rcminputs033101.sav'
help, /structure, DEFINEBF_0603_SW_ASCII
size_result = SIZE( DEFINEBF_0603_SW_ASCII.UT )
imax = size_result[1]  ;=74

print, 'imax',imax
;print,  DEFINEBF_0603_SW_ASCII.UT[0:5]/3600.
imin = 1L
print,  'min',min(DEFINEBF_0603_SW_ASCII.BZ[imin:imax-1]),  ' max',max(DEFINEBF_0603_SW_ASCII.BZ[imin:imax-1])

print, "!Y.CRANGE", !Y.CRANGE, size(!Y.CRANGE)
;1st panel: (2) right; IMF Bz[nT]
scale_factor = (Y_MAX - Y_MIN) / (+50.0 - (-50.0))
offset = (Y_MAX - Y_MIN)*0.50 + Y_MIN
print, 'scale_factor', scale_factor, ' offset', offset
oplot, (DEFINEBF_0603_SW_ASCII.UT[imin:imax-1]/3600. -24.), (DEFINEBF_0603_SW_ASCII.BZ[imin:imax-1]*scale_factor + offset) $
;, yrange=[ Y_MIN, Y_MAX ], ystyle=1  $
;, xrange=[ X_MIN, X_MAX ], xstyle=1  $
;, XTITLE='UT [hours]', YTITLE='IMF Bz [nT]'  $
;, charsize=char_size, charthick=char_thick  ;$
, LINESTYLE=5    ;$

;reference at 0
oplot, (DEFINEBF_0603_SW_ASCII.UT[imin:imax-1]/3600. -24.), (DEFINEBF_0603_SW_ASCII.BZ[imin:imax-1]*0.00 + offset)  $
, LINESTYLE=1    ;$

;plot axis on the right
AXIS, YAXIS=1, YRANGE = (!Y.CRANGE-offset)/scale_factor, YSTYLE = 1, $
   YTITLE = 'IMF B!DZ!N [nT]' $
, charsize=char_size, charthick=char_thick  ;$



; 2nd top panel: SYM-H(Dst)
Y0 = Y00 +(dY +dY_margin)*1.0

print,  'min',min(DEFINEBF_0603_SW_ASCII.DST[imin:imax-1]),  ' max',max(DEFINEBF_0603_SW_ASCII.DST[imin:imax-1])
Y_MIN= -450.00 
Y_MAX=  +50.00 

plot, (DEFINEBF_0603_SW_ASCII.UT[imin:imax-1]/3600. -24.), DEFINEBF_0603_SW_ASCII.DST[imin:imax-1]  $
, yrange=[ Y_MIN, Y_MAX ], ystyle=8  $  ;070306: 1  $
, xrange=[ X_MIN, X_MAX ], xstyle=1  $
;, XTITLE='UT [hours]' $
, YTITLE='SYM-H'  $      ;070306: D!DST!N'  $
, charsize=char_size, charthick=char_thick  $
, LINESTYLE=0    $
;, XTICKINTERVAL=3.0  $
, pos=[X0/X_SIZE, Y0/Y_SIZE, (X0+dX)/X_SIZE, (Y0+dY)/Y_SIZE]  $
,/NOERASE $
, XTICKs=8, xminor=3  $
, XTICKLEN=1, XGRIDSTYLE=1  $ 
, YTICKs=10, Yminor=5 , Ytickname=[' ','-400',' ','-300',' ','-200',' ','-100',' ','0',' ' ]  $
, YTICKLEN=1, YGRIDSTYLE=1 ;$ 

;2nd panel (2) right; P_dyn [n Pa]
print,  'P_dyn: min',min(DEFINEBF_0603_SW_ASCII.P_DYN[imin:imax-1]),  ' max',max(DEFINEBF_0603_SW_ASCII.P_DYN[imin:imax-1])
scale_factor = (Y_MAX - Y_MIN) / (+80.0 - (+0.0))
offset = Y_MIN  ;070306: (Y_MAX - Y_MIN)*0.50 + Y_MIN
print, 'scale_factor', scale_factor, ' offset', offset

oplot, (DEFINEBF_0603_SW_ASCII.UT[imin:imax-1]/3600. -24.), (DEFINEBF_0603_SW_ASCII.P_DYN[imin:imax-1]*scale_factor + offset) $
, LINESTYLE=5 

;plot axis on the right
AXIS, YAXIS=1, YRANGE = (!Y.CRANGE-offset)/scale_factor, YSTYLE = 1, $
   YTITLE = 'P!DDYN!N [n Pa]' $
, charsize=char_size, charthick=char_thick  ;$

;reference
;oplot, (DEFINEBF_0603_SW_ASCII.UT[imin:imax-1]/3600. -24.), (DEFINEBF_0603_SW_ASCII.DST[imin:imax-1]*0.00)  $
;, LINESTYLE=1    ;$


;3rd panel: Jicamarca drift [m/s]
; 3.1. RCM+CTIPe drift
RESTORE, filename="/Users/naomi/idl/plot_ctip/drift/ExBdrift033101_LON283.0.sav"

print, 'sizeVPEQ', sizeVPEQ
time = findgen(sizeVPEQ)*15./60. +12.25 -24.00
print, 'min',min(vpeqX[0:sizeVPEQ-1 ,3]), 'max',max(vpeqX[0:sizeVPEQ-1 ,3])
Y_MAX=  +60.00
Y_MIN=  -Y_MAX


Y0 = Y00 +(dY +dY_margin)*0.0
;QT drift (run1)
plot, time, vpeqX[0:sizeVPEQ-1 ,0]  $
, yrange=[ Y_MIN, Y_MAX ], ystyle=8  $
, xrange=[ X_MIN, X_MAX ], xstyle=1  $
, XTITLE='UT [hours]', YTITLE='Vertical ExB!Cdrift [m/s]'  $
, charsize=char_size, charthick=char_thick  $
, LINESTYLE=3    $
;, XTICKINTERVAL=3.0  $
, pos=[X0/X_SIZE, Y0/Y_SIZE, (X0+dX)/X_SIZE, (Y0+dY)/Y_SIZE]  $
;, color=PlotCol   $
,/NOERASE  $
, XTICKs=8, xminor=3  $
, XTICKLEN=1, XGRIDSTYLE=1  $ 
, YTICKs=6, Yminor=5  $
, YTICKLEN=1, YGRIDSTYLE=1 ;$

;PP drift (run3)
PlotCol= col_min +60.   ;blue;
oplot, time, vpeqX[0:sizeVPEQ-1 ,2]  $
, color=PlotCol   ;$

;DD drift (run2)
oplot, time, vpeqX[0:sizeVPEQ-1 ,1] ;$

;total drift (run4)
PlotCol= col_max -7. ;red
oplot, time, vpeqX[0:sizeVPEQ-1 ,3] $
, color=PlotCol   ;$

;reference at 0
;oplot, time, (vpeqX[0:sizeVPEQ-1 ,3]*0.00)  $
;, LINESTYLE=1    ;$



; 3.2. Dave's drift
restore, filename='/Users/naomi/jastp_special2006/data4JASTP030906/adela_drift033101_slope052306.sav'
size_result = SIZE( JP90_2001.LT_HR ) ;=288
kmax = size_result[1]
print, kmax
;print, 'min',min(JP90_2001.LT_HR[0:kmax-1]), ' max',max(JP90_2001.LT_HR[0:kmax-1])
kmin = 72L  +10L  ;052306:
print, '6LT',JP90_2001.LT_HR[kmin-1:kmin+1]
kkmax = kmax - kmin
print, '18LT',JP90_2001.LT_HR[kkmax-1:kkmax+1]

print, 'min',min(JP90_2001.dH_Drift_East[kmin:kkmax]), ' max',max(JP90_2001.dH_Drift_East[kmin:kkmax])

glonX=283.00 ;jicamarca
ut_hr = findgen(kmax)
ut_hr[0:kmax-1] = JP90_2001.LT_HR[0:kmax-1] - glonX / 15.0 +24.0
PlotCol= col_max-45. ;orange?   -15. ;red
oplot, ut_hr[kmin:kkmax], JP90_2001.dH_Drift_East[kmin:kkmax] $
;, charsize=char_size, charthick=char_thick  ;$
, LINESTYLE=0    $
, thick=4.5      $
, color=PlotCol  ;$

;FS drift
;oplot, ut_hr[kmin:kkmax], JP89_2001.dH_Drift_East[kmin:kkmax]  $
;, LINESTYLE=1    ;$
;oplot, ut_hr[kmin:kkmax], JP91_2001.dH_Drift_East[kmin:kkmax]  $
;, LINESTYLE=3    ;$


Bmag = 2.65E-05 ;[tesla]
offset = 0.00
scale_factor = Bmag * 1.0E+03  ;drift[m/s] --> Efield [mV/m]
;plot axis on the right
AXIS, YAXIS=1, YRANGE = (!Y.CRANGE-offset)*scale_factor, YSTYLE = 1, $
   YTITLE = 'Eastward Electric!CField [mV/m]' $
, charsize=char_size, charthick=char_thick  ;$

DEVICE, /CLOSE
end ;pro ts_plot_fig2_1
