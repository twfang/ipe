;  pro read_EXB $
  pro rd_EXB $
, plot_VEXB_ut, plot_VEXB,mp,lp,NMP,NLP, input_DIR,sw_debug,n_max1

size_results=SIZE(plot_VEXB_ut)
n_max=size_results[1]

VEXB=fltarr(NLP,NMP)
openr, LUN14, input_DIR+'ut_rec', /GET_LUN 
openr, LUN15, input_DIR+'plasma16', /GET_LUN $
, /F77_UNFORMATTED ;$

n=-1L
;while(eof(LUN14) eq 0 ) do begin
while( n lt n_max-1 ) do begin
n=n+1L
 readf, LUN14, record_number, UT_sec
print,n,'rec#',record_number,' UThr',(UT_sec /3600.)
plot_VEXB_ut[n] = UT_sec

mpstop=1L
dum=fltarr(NLP)
for i=0,mpstop-1 do begin
     readu, LUN15, dum ;VEXB
;print, dum
     VEXB[0:NLP-1,i]=dum[0:NLP-1]
endfor

if sw_debug eq 1 then  print, 'VEXB_ms1=',VEXB[lp,mp],mp,lp
for j=0,2 do begin
jth=lp + 10*j
if ( n eq 0 ) then  print ,j,'VEXB: jth=',jth
plot_VEXB[n,j]=VEXB[jth,mp]
endfor ;j
endwhile;(eof(LUN14) eq 0 ) do begin

n_max1=n

FREE_LUN, LUN14
FREE_LUN, LUN15
end;  pro rd_EXB, plot_VEXB_ut, plot_VEXB
