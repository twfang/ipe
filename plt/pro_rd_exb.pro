;  pro read_EXB $
  pro pro_rd_EXB ;$
;, plot_VEXB_ut, plot_VEXB,mp,lp,NMP,NLP, input_DIR,sw_debug, mlat90_2d $
;,n_max1   ;output

input_DIR=$
;'/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r336/trunk/run/ipe_S_29795/'
'/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r319.4/trunk/run/ipe_S_29672/'
sw_exb_zonal=0L
;size_results=SIZE(plot_VEXB_ut)
;n_max=size_results[1]
n_max=1L
NLP=170L
NMP=80L
sw_debug=1L
mlat90=fltarr(NLP)
;mlat90=( -88.1238     -87.6257     -86.2386     -85.5684     -84.3344     -83.5766     -82.4013
;     -81.5840     -80.4295     -79.5634     -78.4094     -77.4983     -76.3318     -75.3764     -74.1877
;     -73.1872     -71.9689     -70.9219     -69.6682     -68.5731     -67.2792     -66.1349     -64.7976
;     -63.6034     -62.2206     -60.9769     -59.5483     -58.2566     -56.7834     -55.4469     -53.9321
;     -52.5554     -51.0042     -49.5940     -48.0135     -46.5780     -44.9773     -43.5264     -41.9163
;     -40.4614     -38.8541     -37.4074     -35.8160     -34.3903     -32.8280     -31.4360     -29.9158
;     -29.6552     -29.3904     -29.1214     -28.8478     -28.5696     -28.2867     -27.9989     -27.7059
;     -27.4076     -27.1038     -26.8539     -26.6002     -26.3424     -26.0805     -25.8143     -25.5437
;     -25.2685     -24.9886     -24.7038     -24.4138     -24.1765     -23.9356     -23.6911     -23.4427
;     -23.1905     -22.9341     -22.6734     -22.4084     -22.1388     -21.8645     -21.6414     -21.4151
;     -21.1853     -20.9521     -20.7152     -20.4746     -20.2300     -19.9814     -19.7286     -19.4714
;     -19.2637     -19.0531     -18.8393     -18.6224     -18.4021     -18.1784     -17.9511     -17.7201
;     -17.4852     -17.2462     -17.0548     -16.8607     -16.6638     -16.4640     -16.2612     -16.0552
;     -15.8460     -15.6333     -15.4171     -15.1973     -14.9049     -14.6057     -14.2990     -13.9845
;     -13.6615     -13.3294     -13.0657     -12.7957     -12.5192     -12.2355     -11.9443     -11.6449
;     -11.4096     -11.1689     -10.9223     -10.6694     -10.4098     -10.1429     -9.93554     -9.72346
;     -9.50625     -9.28361     -9.05508     -8.82022     -8.64002     -8.45573     -8.26707     -8.07380
;     -7.87549     -7.67178     -7.51718     -7.35919     -7.19756     -7.03210     -6.86244     -6.68832
;     -6.55756     -6.42403     -6.28762     -6.14809     -6.00518     -5.85866     -5.74927     -5.63771
;     -5.52389     -5.40754     -5.28861     -5.16683     -4.88908     -4.59414     -4.36233     -4.11728
;     -3.91972     -3.71156     -3.53623     -3.35164     -3.18684     -3.01298     -2.84743     -2.67174
;     -2.49473     -2.30410     -2.10151)

VEXB =fltarr(NLP,NMP)
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
;plot_VEXB_ut[n] = UT_sec

mpstop=NMP
dum=fltarr(NLP)
for i=0,mpstop-1 do begin
print,'mp_read=',i
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
print,mpt,vexbe[0:nlp-1,mpt]
plot, vexbe[0:nlp-1,mpt] $
, yrange=[vmin , vmax ],  ystyle=1 $
, linestyle=0, color=254.
mpt=1L
print,mpt,vexbe[0:nlp-1,mpt]
oplot, vexbe[0:nlp-1,mpt] $
, linestyle=5, color=50.
mpt=80-1L
print,mpt,vexbe[0:nlp-1,mpt]
oplot, vexbe[0:nlp-1,mpt] $
, linestyle=2, color=150.

oplot, mlat90_2d[0,0:nlp-1],vexbe[0:nlp-1,0]*0.0 $
, linestyle=1
stop

;for j=0,2 do begin
;jth=lp; + 10*j
;if ( n eq 0 ) then  print ,j,'VEXB: jth=',jth
;plot_VEXB[n,0]=VEXB[ lp,mp];upward drift
;plot_VEXB[n,1]=VEXBe[lp,mp];zonal drift
;endfor ;j
endwhile;(eof(LUN14) eq 0 ) do begin

n_max1=n

FREE_LUN, LUN14
FREE_LUN, LUN15
if ( sw_exb_zonal eq 1 ) then $
FREE_LUN, LUN17
end;  pro rd_EXB, plot_VEXB_ut, plot_VEXB
