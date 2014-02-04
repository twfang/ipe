;20140108
;purpose: read wind on IPE grid for validation purpose
pro read_wind, ut_hr $
, Vn_ms1 $
;, Un_ms1 $
, luntmp7, luntmp3 $
, MaxFluxTube, nlp, nmp, JMIN_IN,JMAX_IS, sw_debug
;
udum0=fltarr(MaxFluxTube,NLP   )
udum=fltarr(MaxFluxTube,NLP,NMP)
;
for k=1-1,nmp-1 do begin
  if ( k eq 0 ) then begin
    readf, luntmp7, utime_dum
    print, k,'sub-read_wind',utime_dum, utime_dum/3600., ut_hr
  endif
;input to IPE
;readu, luntmp3, Vn_ms1 ;(1:3,1:NPTS,1:nmp) 
;print, 'sub-read_wind',vn_ms1[0,0,0], vn_ms1[0,1116,0]
;output to IPE
  readu, luntmp3, udum0 ;(1:MaxFluxTube,1:NLP,1:nmp,3:3) 
  if ( sw_debug eq 1 ) then print,k, 'sub-read_wind',Udum0[*,130,0]
  for j=1-1,NLP-1 do begin
    for i=1-1,MaxFluxTube-1 do begin
      udum[i,j,k]=udum0[i,j]
    endfor ;i
  endfor ;j
endfor ;k=1-1,nmp-1 do begin

;udum-->vn_ms1
for k=1-1,nmp-1 do begin
   for j=1-1,nlp-1 do begin
      in=JMIN_IN[j]
      is=JMAX_IS[j]

      for i=1-1,is-in+1-1  do begin
         ii = i + (in-1)
         vn_ms1[2-1,ii,k] = udum[i,j,k]
      endfor                    ;i
   endfor                       ;j
endfor                          ;k

;STOP
  if ( sw_debug eq 1 ) then print , 'sub-read_wind finished!'
end;
