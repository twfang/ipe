pro read_fort9001 $
, ph0,ph1,th0,th1

input_DIR='/home/Naomi.Maruyama/p2/ipe/r336.2/trunk/run/ipe_S_7196/'
input_flnm=input_DIR+'fort.9001'
LUN9001=0L
;open
   openr, LUN9001, input_flnm, /GET_LUN
;read
n_read = -1L
mp=0L
lp=0L
while ( eof(LUN9001) eq 0 ) do begin
  n_read = n_read + 1
   readf, LUN9001,mp,lp,phi0,phi1,the0,the1 ;$
;FORMAT=$
;fmt00
;'(A16,2F10.2)'

print,'n_rd=',n_read,' mp=', mp,' lp=',lp,phi0,phi1,the0,the1

ph0[mp-1,lp-1]=phi0
ph1[mp-1,lp-1]=phi1
th0[mp-1,lp-1]=the0
th1[mp-1,lp-1]=the1

endwhile

;close
   FREE_LUN, LUN9001

print, "pro read_fort9001"
end ;pro read_fort9001
