;20140108
;purpose: read wind on IPE grid for validation purpose
pro read_wind, ut_hr $
, Vn_ms1,tn_k,on_m3,n2n_m3,o2n_m3 $
, luntmp7, luntmp3, luntmp8, luntmp9, luntmp10, luntmp11 $
, MaxFluxTube, nlp, nmp, JMIN_IN,JMAX_IS, sw_debug, sw_version_io


sw_version=0L ;wam-ipe neutral output
;sw_version=1L ;debug wind on ipe grid



NMPdum=nmp

;udum0=fltarr(MaxFluxTube,NLP   )
if sw_version eq 0L then begin

if sw_debug eq 1 then  print, 'MaxFluxTube',MaxFluxTube, 'NLP=',NLP,' NMPdum',NMPdum
;   udum =fltarr(MaxFluxTube,NLP,NMPdum,3) ;wind before 20140806
;d   udum =fltarr(MaxFluxTube,NLP,NMPdum,1) ;wind dbg20170421
   udum =fltarr(MaxFluxTube,NLP,NMPdum,3) ;20170824 wind in george's version
   udum1=fltarr(MaxFluxTube,NLP,NMPdum) ;Tn
   udum2=fltarr(MaxFluxTube,NLP,NMPdum) ;o
   udum3=fltarr(MaxFluxTube,NLP,NMPdum) ;n2
   udum4=fltarr(MaxFluxTube,NLP,NMPdum) ;o2

;
;t for k=1-1,nmp-1 do begin
;t   if ( k eq 0 ) then begin
   if sw_version_io eq 0 then begin
      readf, luntmp7, utime_dum
      print, 'sub-read_wind: utime=',utime_dum, utime_dum/3600., ut_hr
      if ( ut_hr ne utime_dum/3600. ) then STOP
   endif                        ;sw_version_io
;t  endif
;input to IPE
;readu, luntmp3, Vn_ms1 ;(1:3,1:NPTS,1:nmp) 
;print, 'sub-read_wind',vn_ms1[0,0,0], vn_ms1[0,1116,0]
;output to IPE

   if sw_version_io eq 0 then begin
      ;readu, luntmp3,  udum      ;wind
      readu, luntmp8,  udum1     ;tn
      readu, luntmp9,  udum2     ;on_m3
      ;readu, luntmp10, udum3       ;n2n_m3
      ;readu, luntmp11, udum4       ;o2n_m3
   endif else if sw_version_io eq 1 then begin
      readu, luntmp3,  udum1      ;tn
      readu, luntmp3,  udum       ;un_ms1
      readu, luntmp3,  udum2      ;on_m3
      readu, luntmp3,  udum3      ;n2n_m3
      readu, luntmp3,  udum4      ;o2n_m3
   endif ;sw_version_io



;dbg20140806
;   print,'check NORTHward wind', MAX(udum[*,*,*,1]), MIN(udum[*,*,*,1])
   print,'check Parallel wind', MAX(udum[*,*,*,0]), MIN(udum[*,*,*,0])
   print,'read check Tn', MAX(udum1), MIN(udum1)
   print,'check [O]', MAX(udum2), MIN(udum2)

;udum-->vn_ms1
   for k=1-1,nmp-1 do begin
      for j=1-1,nlp-1 do begin
         in=JMIN_IN[j]
         is=JMAX_IS[j]
         
         for i=1-1, is-in+1-1  do begin
            ii = i + (in-1)
;dbg if k eq 0 AND i ge 0 AND i le 5 then  print, i,ii,j,in,k
            
;nm20141015: positive eastward
           vn_ms1[1-1,ii,k] = udum[i,j,k,1-1]
;nm20141015: positive northward
           vn_ms1[2-1,ii,k] = udum[i,j,k,2-1]
;positive upward 
           vn_ms1[3-1,ii,k] = udum[i,j,k,3-1]
;nm20141015: positive UPward
;t           vn_ms1[3-1,ii,k] = udum[i,j,k,3-1]
           tn_k[ii,k]   = udum1[i,j,k]
           on_m3[ii,k]  = udum2[i,j,k]
           n2n_m3[ii,k] = udum3[i,j,k]
           o2n_m3[ii,k] = udum4[i,j,k]
            


         endfor                 ;i
      endfor                    ;j
   endfor                       ;k
   

print,'check on_m3', MAX( on_m3 ), MIN( on_m3 )

mp=0L
print, '(1)udum2',udum2[0:4,0,mp]
print, '(1)on_m3',on_m3[0:4  ,mp]

print, '(2)udum2',udum2[   0:4   ,1,mp]
print, '(2)on_m3',on_m3[1115:1119  ,mp]


mp=79L
print, '(3)udum2',udum2[0:4,0,mp]
print, '(3)on_m3',on_m3[0:4  ,mp]

print, '(4)udum2',udum2[   0:4   ,1,mp]
print, '(4)on_m3',on_m3[1115:1119  ,mp]
;STOP
endif else if sw_version eq 1L then begin
;t print,'(A)'
   udum1=fltarr(MaxFluxTube,NLP,NMP,3) ;runs after 20140806
;print,'(B)'   
;t   readu, luntmp3, udum         ;(1:MaxFluxTube,1:NLP,1:nmp,3:3) ;un_ms1
   readf, luntmp3, udum1 ;ascii
;d print,'(C)'

;udum-->vn_ms1
   for k=1-1,nmp-1 do begin
      for j=1-1,nlp-1 do begin

;d print,'(D)',k,j
         in=JMIN_IN[j]-1
         is=JMAX_IS[j]-1
         
;d print,'(E)',in,is
         for i=1-1, is-in+1-1  do begin
            ii = i + in
            
;d print,'(F)',i,ii

            vn_ms1[1-1,ii,k] = udum1[i,j,k,1-1] ;positive eastward
            vn_ms1[2-1,ii,k] = udum1[i,j,k,2-1] ;positive northward???
         endfor                 ;i
      endfor                    ;j
   endfor                       ;k

endif ;sw_version eq 0L then begin

  if ( sw_debug eq 1 ) then print , 'sub-read_wind finished!'
end;
