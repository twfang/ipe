;  pro read_EXB $
  pro rd_EXB $
, plot_VEXB_ut, plot_VEXB,mp,lp,NMP,NLP, input_DIR,sw_debug, mlat90_2d $
,n_max1   ;output

;input_DIR='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r336/trunk/run/ipe_S_29795/'
sw_exb_zonal=0L
size_results=SIZE(plot_VEXB_ut)
n_max=size_results[1]
;n_max=1
;NLP=170L
;NMP=80L
;sw_debug=1L

VEXB=fltarr(NLP,NMP)
VEXBe=fltarr(NLP,NMP);zonal drift
openr, LUN14, input_DIR+'ut_rec', /GET_LUN 
openr, LUN15, input_DIR+'plasma16', /GET_LUN $
, /F77_UNFORMATTED
if ( sw_exb_zonal eq 1 ) then $
openr, LUN17, input_DIR+'plasma17', /GET_LUN $
, /F77_UNFORMATTED

n=-1L
;while(eof(LUN14) eq 0 ) do begin
while( n lt n_max-1 ) do begin
n=n+1L
 readf, LUN14, record_number, UT_sec
print,n,'rec#',record_number,' UThr',(UT_sec /3600.),' UTs',UT_sec
plot_VEXB_ut[n] = UT_sec

mpstop=NMP
dum=fltarr(NLP)
for i=0,mpstop-1 do begin
;print,'mp_read=',i
   readu, LUN15, dum            ;VEXB upward
;print, dum
   VEXB[ 0:NLP-1,i]=dum[0:NLP-1]
   if ( sw_exb_zonal eq 1 ) then begin
      readu, LUN17, dum         ;VEBeast
      VEXBe[0:NLP-1,i]=dum[0:NLP-1]
   endif
endfor

if sw_debug eq 1 then  begin
   print, 'VEXB_ms1: upward',VEXB[lp,mp],mp,lp,' eastward=',VEXBe[lp,mp]
   if ( sw_exb_zonal eq 1 ) then print, 'VEXB_ms1: eastward=',VEXBe[lp,mp]
endif

;20140225
iwindow=2L
DEVICE, RETAIN=2, DECOMPOSED=0
WINDOW,iwindow,XSIZE=800,YSIZE=800
!p.multi=[0,1,1,0]
vmax=200.
vmin=vmax *(-1.)
mpt=0L
;print,mpt,vexbe[0:nlp-1,mpt]
plot, vexbe[0:nlp-1,mpt] $
, yrange=[vmin , vmax ],  ystyle=1 $
, linestyle=0, color=254.
mpt=1L
;print,mpt,vexbe[0:nlp-1,mpt]
oplot, vexbe[0:nlp-1,mpt] $
, linestyle=5, color=50.
mpt=80-1L
;print,mpt,vexbe[0:nlp-1,mpt]
oplot, vexbe[0:nlp-1,mpt] $
, linestyle=2, color=150.

oplot, mlat90_2d[0,0:nlp-1],vexbe[0:nlp-1,0]*0.0 $
, linestyle=1
stop

;for j=0,2 do begin
;jth=lp; + 10*j
;if ( n eq 0 ) then  print ,j,'VEXB: jth=',jth
plot_VEXB[n,0]=VEXB[ lp,mp];upward drift
plot_VEXB[n,1]=VEXBe[lp,mp];zonal drift
;endfor ;j
endwhile;(eof(LUN14) eq 0 ) do begin

n_max1=n

FREE_LUN, LUN14
FREE_LUN, LUN15
if ( sw_exb_zonal eq 1 ) then $
FREE_LUN, LUN17
end;  pro rd_EXB, plot_VEXB_ut, plot_VEXB
