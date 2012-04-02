  pro read_EXB, plot_VEXB_ut, plot_VEXB,mp,lp,NMP,NLP, input_DIR,sw_debug,n_max1

size_results=SIZE(plot_VEXB_ut)
n_max=size_results[1]

VEXB=fltarr(NMP,NLP)
openr, LUN14, input_DIR+'ut_rec.log', /GET_LUN 
openr, LUN15, input_DIR+'plasma16', /GET_LUN $
, /F77_UNFORMATTED ;$

n=-1L
while(eof(LUN14) eq 0 ) do begin
n=n+1L
 readf, LUN14, record_number, UT_sec
print,n,'record_number=',record_number,' UT_hr=',(UT_sec /3600.)
plot_VEXB_ut[n] = UT_sec

     readu, LUN15, VEXB
if sw_debug eq 1 then  print, 'VEXB_ms1=',VEXB[mp,lp],mp,lp
plot_VEXB[n]=VEXB[mp,lp]
endwhile;(eof(LUN14) eq 0 ) do begin

n_max1=n

FREE_LUN, LUN14
FREE_LUN, LUN15
end;  pro read_EXB, plot_VEXB_ut, plot_VEXB
