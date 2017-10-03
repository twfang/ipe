;;;i must include this routine into rc_input_txt.pro
pro rd_cpcp

runid='2015' ;2013';quiet';
sw_output2png=1L
sw_save=0L
sw_cmp=1L

if runid eq '2013' then begin
  rundir1='1461312397'
  rundir2='1461343191'
  path='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/mpi20160330v2/run'
  nLevPI1=1078L ;00-17ut; <583020
  nLevPI2=371L  ;17-24ut:  583080<utime<604800
endif else if runid eq '2015' then begin
  rundir1='1461417243'
  rundir2='1461436195'
  path='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/mpi20160330v2/run1'
  nLevPI1=1020L ;00-17ut; <579540
  nLevPI2=421L  ;17-24ut:  579600<utime<604800
endif else if runid eq 'quiet' then begin
  rundir1='1462618349'
  rundir2='1462657622'
  path='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/mpi20160330v2/run2'
  nLevPI1=1020L ;00-17ut; <579540
  nLevPI2=421L  ;17-24ut:  579600<utime<604800
endif
rundir0='_ipe_theia_intel_parallel2_93'
nLevPI=1440L 
flnm='fort.1002'
flnm_png='cpcp'+runid+'0317.png'
lun0=0L
cpcp=fltarr(nLevPI)
uthr=fltarr(nLevPI)
ut0=518400L
ut=0L

;file1 00--17ut
openr, LUN0, path+'/'+rundir1+rundir0+'/'+flnm, /GET_LUN
for i=0,nlevpi1-1 do begin
   readf, LUN0, ut,potmin,potmax,dum
   print,i   , ut,potmin,potmax,dum
   cpcp[i]=dum
   uthr[i]=(ut-ut0)/3600.
endfor
free_lun, lun0

;file2 17--24ut
openr, LUN0, path+'/'+rundir2+rundir0+'/'+flnm, /GET_LUN
for i=nlevpi1, (nlevpi1+nlevpi2-1) do begin
   if i ge nLevPI then break
   readf, LUN0, ut,potmin,potmax,dum
   print,i   , ut,potmin,potmax,dum
   cpcp[i]=dum
   uthr[i]=(ut-ut0)/3600.
endfor
free_lun, lun0
;print,'cpcp=',cpcp

iwindow=1L
DEVICE, RETAIN=2, DECOMPOSED=0
WINDOW,iwindow,XSIZE=1200,YSIZE=1000

loadct,0

x_min=0.
x_max=24.
y_max=+280.;210. ;22.
y_min=10.;-y_max
char_size=2.0
char_thick=2.0

if sw_cmp eq 0 then $
  titleMain='PCP: '+runid+'-03-17' $
else if sw_cmp eq 1 then $
  titleMain='PCP Comparison: '+runid+'-03-17  vs. 2013-03-17'


print,'MAX PCP=',max(cpcp,max_subscript),uthr(max_subscript)
print,'MIN PCP=',min(cpcp,min_subscript),uthr(min_subscript),(uthr(min_subscript)*3600.)
plot,uthr,cpcp $
  ,xrange=[x_min,x_max], xstyle=1  $
  ,yrange=[y_min,y_max], ystyle=1  $
  ,title= titleMain $
  ,xtitle='UT [hr]' $
  ,ytitle='PCP [kV]' $
  , charsize=char_size, charthick=char_thick $
  ,/NODATA

loadct,39
red=230
green=150
oplot,uthr,cpcp $
,color=red

if sw_cmp eq 1 then begin
   restore,'cpcp2013.sav'
   oplot,uthr,cpcp $
         ,color=green

;quiet time : cpcp=44.14kV
   restore,'cpcpquiet0317.sav'
loadct,0
   oplot,uthr,cpcp $
         , linestyle=1

loadct,39
   x1=0.9
   xyouts,x1,0.8,'2015',color=red $
          ,charsize=char_size,charthick=char_thick $
          , /norm,/noclip
   xyouts,x1,0.35,'2013',color=green $
          ,charsize=char_size,charthick=char_thick $
          , /norm,/noclip

loadct,0
   xyouts,x1,0.20,'quiet' $ ;,color=green $
          ,charsize=char_size,charthick=char_thick $
          , /norm,/noclip

   flnm_png='comp_cpcp'+flnm_png
  ;output_png,'comp'+flnm_png
endif ; sw_cmp eq 1 then begin


if sw_output2png eq 1 then  output_png,flnm_png

if sw_save eq 1 then  save, filename='cpcp'+runid+'0317.sav'
print,'pro rd_cpcp finished successfully'

end;pro rd_cpcp
