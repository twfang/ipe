pro read_sunlons $
, sunlons1,TEST,rundir, LUN2013, n_read,n_read_freq

  if n_read eq 0 then begin
;  input_DIR='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/runs/'+TEST+'/trunk/run/'+rundir
;  input_DIR='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/'+TEST+'/run/'+rundir
;  input_DIR='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/'+TEST+'/run1/'+rundir
  input_DIR='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/'+TEST+'/run2/'+rundir
  print, 'read_fort2013:input_DIR=', input_DIR
  input_flnm=input_DIR+'/fort.2013'
;LUN2013=0L
;i need to think about what to do after second time...
;open

   openr, LUN2013, input_flnm, /GET_LUN
endif

;internal read loop
inRead = 0L
inReadMax = n_read_freq

while ( inRead lt inReadMax ) do begin
   inRead = inRead + 1
   readf, LUN2013,sunlons1

print,'internal n_read sunlon=',inRead,' sunlons1=', sunlons1

endwhile

;i need to think about what to do after second time...
;close
;   FREE_LUN, LUN2013

print, "pro read_sunlons"
end ;pro read_sunlons
