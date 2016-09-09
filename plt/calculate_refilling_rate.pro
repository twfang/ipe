pro calculate_refilling_rate, den, ut_hr_save, n_read_max,rate,sw_debug

i0=1
i1=n_read_max-1
t0 = ut_hr_save[i0] ;hrs
t1 = ut_hr_save[i1] ;hrs
dt = (t1 - t0) ;=24hrs
if sw_debug eq 1 then   print,'t0=', t0,' t1=', t1,' dt=', dt

n0 = den[i0] ;cm-3
n1 = den[i1] ;cm-3
;dn = n1 - n0
if sw_debug eq 1 then   print,'n0=', n0,' n1=',n1,' dn=',(n1-n0)

;refilling rate
if ( dt gt 0.00000 ) then begin
   rate = (n1-n0) / ( dt/24. )
if sw_debug eq 1 then    print, 'rate [cm-3 day-1] =', rate
endif else  begin
   print, 'STOP! INVALID dt=', dt
   STOP
endelse

end ;pro calculate_refilling_rate
