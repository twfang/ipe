;20140108
;purpose: read wind on IPE grid for validation purpose
pro read_wind, ut_hr, Vn_ms1, luntmp7, luntmp3

readf, luntmp7, utime_dum
print, 'sub-read_wind',utime_dum, utime_dum/3600., ut_hr

readu, luntmp3, Vn_ms1 ;(1:3,1:NPTS,1:nmp) 
print, 'sub-read_wind',vn_ms1[0,0,0], vn_ms1[0,1116,0]

;STOP
end;
