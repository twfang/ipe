pro dbg_frict_heating

sw_plot_Tei=0L; 1
sw_frictHeat=1L
;titlemplp='Mp4Lp52'
;titlemplp='Mp4Lp51'
titlemplp='Mp3Lp52'
if sw_frictHeat eq 1 then begin
  titleMain='sw_frictHeat1'
  ;rundir0='1492015228'
  ;rundir0='1492802085'
  rundir0='1493047632'
endif else if sw_frictHeat eq 0 then begin
  titleMain='sw_frictHeat0'
;  rundir0='1492015341'
;  rundir0='1492802281'
  rundir0='1493047761'
endif
rundir1=rundir0+'_ipe_theia_intel_parallel2_93'
path='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/mpi20160330v4/run/'+rundir1+'/'

openr, lun, path+'fort.8011', /get_lun

tmp_string=' Ti: ALT      Te       Ti       Tn      F        dT/dt      C_ie      C_in    HFLUX_UP  HFLUX_LO   DIVV_UP  DIVV_LO    DIVCOND  DIVCONV             FHEAT     FRICT1     FRICT2     FRICT3'

readf, lun, tmp_string
print, tmp_string

;tmp_string2='mp   4lp  52j   6'
tmp_string2='mp   4lp  51j   6'


jmax=105L
alt=fltarr(jmax)
Te=fltarr(jmax)
Ti=fltarr(jmax)
Tn=fltarr(jmax)
Ftotal=fltarr(jmax)
dTi_over_dt=fltarr(jmax)
C_ie=fltarr(jmax)
C_in=fltarr(jmax)
Hflux_up=fltarr(jmax)
Hflux_lo=fltarr(jmax)
Fheat=fltarr(jmax)

ncnt=1L
;j loop
while ( ncnt le jmax) do begin


   readf, lun, tmp_string2
   print, tmp_string2


;ALT  Te   Ti   Tn   F    dT/dt  C_ie   C_in   HFLUX_UP HFLUX_LO DIVV_UP DIVV_LO  DIVCOND  DIVCONV FHEAT     FRICT1     FRICT2     FRICT3'
;Z(J),TI(3)TI(2)TN(J)F(1) DERIVT ElossI,IlossN,HFLUX_UP,HFLUX_LO,DIVV_UP,DIVV_LO, HFLUXU-L,DIVVU-L,(HFLUX_UP-HFLUX_LO)/(1.0E-33+DIVV_UP-DIVV_LO),  fheat,fheat/F(1)
   readf, lun, $
dum0, dum1,dum2,dum3,dum4,dum5,  dum6,  dum7,  dum8,    dum9    ,dum10,  dum11,   dum12,   dum13,   dum14,                                         dum15,dum16 $
;,FORMAT='(4F9.1,1P,22E10.2)'
          ,FORMAT='(4F9.1,13E10.2)'

alt[ncnt-1]=dum0
Te[ncnt-1] =dum1
Ti[ncnt-1] =dum2
Tn[ncnt-1] =dum3
Ftotal[ncnt-1] =dum4
dTi_over_dt[ncnt-1] =dum5
C_ie[ncnt-1] =dum6
C_in[ncnt-1] =dum7
Hflux_up[ncnt-1] =dum8
Hflux_lo[ncnt-1] =dum9
Fheat[ncnt-1] =dum15

   print,ncnt, dum0,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9 $
         ,dum10,dum11,dum12,dum13,dum14,dum15,dum16 $
         ,FORMAT='(i5,4F9.1,13E10.2)'

   ncnt=ncnt+1

endwhile
; end j loop
print, ncnt

free_lun, lun

;plot

  device, decomposed = 0 ,retain=2
  window, 0 ,XSIZE=1000,YSIZE=1000
  n_ldct=39
  loadct, n_ldct

;!P.MULTI=[0,2,1,0,1]
;1:# of plot columns
;2:# of rows


x_min=200.
x_max=1300.

y_min=100.
y_max=1300.

x_min1= -0.000281000; -1.0E-3
x_max1= +1.0E-3


if sw_plot_Tei eq 1 then begin
plot,  Ti,alt $
  ,xrange=[x_min,x_max], xstyle=1  $
  ,yrange=[y_min,y_max], ystyle=1  $
  ,title=titleMain+': '+rundir1 $
  ,linestyle=0
oplot, Te,alt,linestyle=2
oplot, Tn,alt,linestyle=3
endif ;sw_plot_Tei eq 1 then begin

;dbg20170421
print,titleMain
for i=0,ncnt-2  do begin
   if alt[i] ge 120. AND alt[i] lt 450. then begin
      rat=fheat[i]/ftotal[i]
      print, i,alt[i],Ftotal[i],Fheat[i],rat, format='(" i=",i4,"  alt=",f5.0,"  Ftotal=",e11.2,"  Fheat=",e11.2, f8.2)'
   endif
endfor
print,'ftotal MIN=', MIN(ftotal), MAX(ftotal)
;
plot,  Ftotal,alt $
 ; ,xrange=[x_min1,x_max1], xstyle=1, /xlog  $
  ,xrange=[-1.0E-4,x_max1], xstyle=1, /xlog  $
  ,yrange=[y_min,y_max], ystyle=1  $
  ,title=titleMain+': '+rundir1+tmp_string2 $
  ,linestyle=0, thick=5
;,psym=4

oplot, dTi_over_dt,alt,linestyle=1
oplot, C_ie,alt,linestyle=2, color=60., thick=2
oplot, C_in,alt,linestyle=3, color=160., thick=2
oplot, Hflux_up,alt,linestyle=4
oplot, Hflux_lo,alt,linestyle=5
print,'fheat MIN=', MIN(fheat), MAX(fheat)
oplot, Fheat,alt,linestyle=6, thick=5, color=250.

;ref
alt=alt*0.+292.
fheat=findgen(jmax)*.5  + 1.0E-15
oplot, Fheat,alt,linestyle=1, thick=1

if sw_plot_Tei eq 0 then $
  flnm_png=titleMain+'_'+rundir0+'_'+titlemplp+'_Terms.png' $
else if sw_plot_Tei eq 1 then $
  flnm_png=titleMain+'_'+rundir0+'_'+titlemplp+'.png'

output_png,flnm_png
print,'pro dbg_frict_heating finished successfully'
end ;pro dbg_frict_heating
