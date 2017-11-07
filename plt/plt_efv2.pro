;20140221: sw_180=1 gives strange plots for the first 4 plots. why
;only the first 4? why not all of them?
;20121120: version 3: ed190 dimension modified to (2*nlp,nmp)
;20120530: fort.2010 added for ed2
;pro plot_efield ;shorter name
pro plt_efv2
;runDATE='20170727'
;runDATE='20170917'
runDATE='20171025'
sw_debug=0L
sw_output2file=1L ;1'PNG' ;0NONE';
sw_output2save=1L ;1save ;0NONE';
sw_plt_cntr=1L ;1:rectangular; 2:polar
sw_plt_exb=0L ;1
version=3  ;2: 20120530; 3:20121120
sw_180=0L ;1:-180<+180; 0:0~360
title_res=$
'low20120709'
;'2xdyn';'
;low';dyn';'low' ; 'high'
utime_min=0;864000;550800.;547200.;###CHANGE
utime_max=utime_min;+600.;to plot ExB time variation on the 6th panel 
;t if sw_plt_exb eq 0 then begin 
  iplot_max=6-1L 
;t  freq_plot_sec=900
;t endif else if sw_plt_exb eq 1 then begin 
;t  iplot_max=6-2L
  freq_plot_sec=utime_max-utime_min
;t endif
;runDIR='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/runs/';before
;runDIR='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/' ;after
runDIR0='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/' ;20170916
figDIR='/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/fig/efield/'

;TEST1='JQIPEr420' ;ref20170727
;TEST1='IPEOptimization' ;testing20170917
;TEST1='ed1130Issue' ;testing20170918
TEST1='raw_high_lat' ;debuging20171025

parallelism='parallel_1'
runid='1508957320'
;parallelism='serial'
;runid='1508956570'

TEST0='20171025IpeOptimization/ipe/run/'
;TEST2='debugE20160401/run/1459548926_ipe_theia_intel_serial2' ;after
;TEST2='tmp20151117/trunk/run4/1459553660_ipe_theia_intel_serial2' ;before
;TEST2='tmp20151117/trunk/run4/1459554325_ipe_theia_intel_serial2' ;after
;TEST2='20170726mergeIPEOptimization2SwpcIpeTest/JQIPEr420/IPEr420/run/1501180820_ipe_theia_intel_parallel_464' ;20170727
;TEST2='20170811testDevelopmentV3/ipe/run/1505642824_ipe_theia_intel_parallel_464' ;20170916
;TEST2='20170811testDevelopmentV3/ipe/run/1505719032_ipe_theia_intel_serial' ;20170918
;TEST2='20170811testDevelopmentV3/ipe/run/1505733014_ipe_theia_intel_parallel_464' ;20170918;module_interp_poten.f
TEST2=TEST0+runid+'_ipe_theia_intel_'+parallelism ;20171025 raw_high_lat test

plot_DIR=$
;'/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/fig/ef/20150317/after20160401/'
;'/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/fig/ef/20150317/before20151210/'
;'/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/fig/20170916/'+TEST1+'/' ;20170916
figDIR+runDATE+'/'+TEST1+'/' ;20171025
print,' plot_DIR=', plot_DIR
input_DIR=runDIR0+TEST2+'/'


utime=0L
freq_output=5 ;min
n_read_max=3L; ###CHANGE
utsec_save=fltarr(n_read_max)
utsec_save[*]=-9999L


nmlon=180L
nmlat=90L
poten=fltarr(nmlon+1,nmlat+1)
ed1130=fltarr(nmlon+1,nmlat+1)
ed2130=fltarr(nmlon+1,nmlat+1)
ylatm=fltarr(nmlat+1)
mlat130=fltarr(nmlat+1)
mlon130=fltarr(nmlon+1)
mlat90_0=fltarr(nmlat+1)
nmp=80L
;nmp=1L  ;dbg20141111
if ( title_res eq 'low' ) OR ( title_res eq 'low20120709' ) then $
   nlp=170L $
else if ( title_res eq '2xdyn' ) then $
   nlp=93L $
else if ( title_res eq 'dyn' ) then $
   nlp=45L                      ;dyn

mlat90_1=fltarr(nlp*2)
ed190=fltarr(nmp,nlp*2) ;mV/m
ed190_save=fltarr(n_read_max) 

ed290=fltarr(nmp,nlp*2)
ed290_save=fltarr(n_read_max)

mlon90_2d=fltarr(nmp,nlp*2)
mlat90_2d=fltarr(nmp,nlp*2)
formatE='(20E12.4)'
formatF='(20f10.4)'
formatF1='(i4,f10.4)'
formatI='(i12)'

openr, LUN10, input_DIR+'fort.2010', /GET_LUN ;utime
;openr, LUN10, input_DIR+'fort.2009', /GET_LUN ;20120125utime
openr, LUN0, input_DIR+'fort.2000', /GET_LUN ;potent
openr, LUN1, input_DIR+'fort.2001', /GET_LUN ;ed1
openr, LUN2, input_DIR+'fort.2002', /GET_LUN ;ed2
openr, LUN3, input_DIR+'fort.2003', /GET_LUN ;ylatm
readf, LUN3, ylatm,  FORMAT=formatF
mlat130=ylatm-90.
openr, LUN4, input_DIR+'fort.2004', /GET_LUN ;ylonm
readf, LUN4, mlon130,  FORMAT=formatF

if ( sw_180 eq 1L ) then begin
for j=0,nmlon+1-1 do begin
  if ( mlon130[j] ge 180. ) then  mlon130[j]=mlon130[j]-360.
endfor
endif

openr, LUN7, input_DIR+'fort.2007', /GET_LUN ;theta90
readf, LUN7, mlat90_0,  FORMAT=formatF
openr, LUN8, input_DIR+'fort.2008', /GET_LUN ;ed1_90
openr, LUN9, input_DIR+'fort.2009', /GET_LUN ;ed2_90
;20130408: fort.2006 does not exist any more!!!
get_gl, mlat90_1,title_res,sw_debug
;openr, LUN6, input_DIR+'fort.2006', /GET_LUN
;string='GL'
;readf, LUN6 , string
;print, string
;for i=1,nlp do begin
;readf, LUN6, ii,lat,  FORMAT=formatF1
;mlat90_1(ii-1)=lat ;NH
;readf, LUN6, ii,lat,  FORMAT=formatF1
;mlat90_1(ii-1)=lat ;SH
;endfor

if sw_debug eq 1L then begin
  print,'after mlat90_1'
 for lp=0,nlp do  print,(lp+1),mlat90_1(lp)
;stop
endif
dlonm90km=4.50 ;deg
mlon90=findgen(nmp)*dlonm90km
if sw_debug eq 1 then $
   print,'mlon90',mlon90
  for i=0,nlp*2-1 do begin
     for mp=0,nmp-1 do begin
        mlon90_2d[mp,i]=mlon90[mp]

        if ( sw_180 eq 1L ) then begin
           if ( mlon90_2d[mp,i] ge 180. ) then mlon90_2d[mp,i]=mlon90_2d[mp,i]-360. 
        endif


      endfor;mp
   endfor   ;i

for j=0,nmp-1 do begin
mlat90_2d[j,0:nlp*2-1]=mlat90_1[0:nlp*2-1]
endfor
;dbg20140225
;print,'mlat90_1=',mlat90_1[0:nlp]
;STOP

n_read=-1L
while(eof(LUN10) eq 0 ) do begin
n_read=n_read+1
readf, LUN10, utime,  FORMAT=formatI
print,'n_read=', n_read,' utime=', utime
utsec_save[n_read]=utime
print,'rd#',n_read,' uts',utime ;,' uts_save',utsec_save[n_read]

readf, LUN0, poten,  FORMAT=formatE
readf, LUN1, ed1130,  FORMAT=formatE
readf, LUN2, ed2130,  FORMAT=formatE
if ( version eq 2 ) then begin  ;nm201205
  readf, LUN8, ed190,  FORMAT=formatE  
  readf, LUN9, ed290,  FORMAT=formatE
endif else $
if ( version eq 3 ) then begin  ;nm20121120
  edum=fltarr(2,nlp,nmp)
  readf, LUN8, edum,  FORMAT=formatE  ;ed190
  for mp=0,nmp-1 do begin
    ii=-1
    for lp=0,nlp-1 do begin
      ii=ii+1
      ed190[mp,ii] = edum[1,lp,mp] ;SH
    endfor ;lp
    for lp=nlp-1,0,-1 do begin
      ii=ii+1
      ed190[mp,ii] = edum[0,lp,mp];NH
    endfor ;lp
  endfor ;mp

  readf, LUN9, edum,  FORMAT=formatE   ;ed290
  for mp=0,nmp-1 do begin
    ii=-1
    for lp=0,nlp-1 do begin
      ii=ii+1
      ed290[mp,ii] = edum[1,lp,mp];SH
    endfor ;lp
    for lp=nlp-1,0,-1 do begin
      ii=ii+1
      ed290[mp,ii] = edum[0,lp,mp];NH
    endfor ;lp
  endfor ;mp
endif ;( version eq 3 ) then begin  ;nm20121120
mp=0L
;mp=17L ;st santin
;mp=4-1L
;mp=8-1L
;mp=60-1L
if ( title_res eq 'low' ) OR ( title_res eq 'low20120709' )  then $
  lp=130-1L $  ;-9.05[deg] mp=0;jicamarca: ht=254.107km
else if ( title_res eq '2xdyn' ) then $
  lp=66-1L  $ ;-9.506[deg] mp=0;jicamarca: ht=???km
else if ( title_res eq 'dyn' ) then $
  lp=34-1L  ;dyn-9.05[deg] mp=0;jicamarca: ht=254.107km

;lp=45L  ;-31.43605[deg];mp=0  Arecibo:31; ht=2504.20km
;lp=28L  ;-56.78335[deg];mp=0 MSH:57; ht=15159.6km
;lp=39L  ;-40.4614[deg]; st santin:40; ht=4790.30km
;lp=46L  ;-29[deg]

ed190_save[n_read]=ed190[mp,lp]
ed290_save[n_read]=ed290[mp,lp]
if ( n_read eq 0 ) then  print,'plt_ef: mp=',mp,' lp=',lp,' mlat90=',mlat90_2d[mp,lp]
title_exb='lat='+STRTRIM( string(mlat90_2d[mp,lp], FORMAT='(F6.2)'),1 )+' mp='+STRTRIM( string((mp+1), FORMAT='(i2)'),1 ) ;IDL convention

if ( utime gt utime_max ) then BREAK ;exit from while read loop
;if $
; ( utime gt utsec_save[0]) AND $
;;( utime ge 605700 ) and $  ;tmp20121120
; ( ( (utime-utsec_save[0]) MOD freq_plot_sec ) LT 0.00001 ) then begin




;20140225 separated out from plt_efv2.pro
if ( sw_plt_cntr eq 1 ) then $
   plt_cntr_fill $
 , iplot_max,mlon90_2d,mlat90_2d, sw_180,mlat130,poten,ed1130,ed2130,ed190,ed290,sw_debug,mlon130,mlat90_0,utime,runDATE,TEST2,plot_DIR,mp,lp, sw_output2file,sw_output2save, parallelism $

else if ( sw_plt_cntr eq 2 ) then $
   plt_cntr_fill_plr $
 , iplot_max,mlon90_2d,mlat90_2d, sw_180,mlat130,poten,ed1130,ed2130,ed190,ed290,sw_debug,mlon130,mlat90_0,utime,runDATE,TEST2,plot_DIR,mp,lp, sw_output2file


;plot ExB time variation
 if ( sw_plt_exb eq 1 ) and ( utime eq utime_max ) then begin
;20140225 separated out
   plt_drft $
 ,mp,lp,NMP,NLP, input_DIR,sw_debug,n_read , utsec_save,ed190_save, title_exb,plot_DIR,utime,runDATE,mlat90_2d,sw_output2file,utime_min
 endif                           ;if ( sw_plt_exb eq 1 ) and ( utime eq utime_max ) then begin???


;endif ;( ( (utime-utime_save) MOD freq_plot_sec ) LT 0.00001 ) then begin

endwhile ;(eof(LUN10) eqq 0 ) do begin


free_lun,lun10
free_lun,lun0
free_lun,lun1
free_lun,lun2
free_lun,lun3
free_lun,lun4
;free_lun,lun6
free_lun,lun7
free_lun,lun8
free_lun,lun9

;print,'mlat90='
;for i=0,nlp*2-1 do print,i,mlat90_2d[0,i]

end ;pro plt_efv2
