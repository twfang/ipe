pro read_plasma_bin,LUN,UT_hr, XIONN_m3,XIONV_ms1,TE_TI_k,VEXB,sw_debug $
,sw_3DJ,je_3d,sw_hr,hrate, sw_dif, sw_lun $
,NMP_in, sunlons1 $
, sza_rad



UT_sec=0L
record_number=0L
; read plasma1
if ( sw_lun[0] eq 1 ) then begin
     readf, LUN[0], record_number, UT_sec
 UT_hr = UT_sec /3600.
print,' rec#',record_number,' UThr=',UT_hr,' UTs=', UT_sec
endif ;( sw_lun[0] eq 1 ) then begin

size_result=size(XIONN_m3)
if ( sw_debug eq 1 ) then  $
   print,'read_plasma: size=',size_result

NPTS2D=size_result[2]
NMP=size_result[3]
if ( sw_debug eq 1 ) then  $
  print, 'read_plasma: NMP',NMP
if ( NMP_in eq 1 ) then NMP=NMP_in
dum=fltarr(NPTS2D,NMP)

;
size_result=size(VEXB)
NLP=size_result[2]
if ( sw_debug eq 1 ) then  $
  print, 'read_plasma:NLP=',NLP

if ( sw_lun[2] eq 1 ) or (sw_lun[20] eq 1)  then begin
   jth=0                        ;o+

   if ( sw_lun[2] eq 1 ) then $
      readu, LUN[2], dum $
   else if ( sw_lun[20] eq 1 ) then $ 
      readu, LUN[20], dum
 



XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
if sw_debug eq 1 then  begin
  print, 'read_plasma:o+ XIONN_m3=',jth,XIONN_m3[jth,60,0]

;dbg20151116 detect nan
  indicesN = where( dum eq !values.f_nan, countN )
  indicesZ = where( dum eq 0.0, countZ )
  if countN ne 0 then begin
     print, 'o+ Nan',countN,indicesN
     STOP
  endif
  if countZ ne 0 then begin
     print, 'o+ Zero',countZ,indicesZ
     STOP
  endif
endif ;sw_debug eq 1 then  begin
endif ;( sw_lun

if ( sw_lun[6] eq 1 ) or (sw_lun[20] eq 1) then begin
   jth=1                        ;h+

  if ( sw_lun[6] eq 1 ) then $
     readu, LUN[6], dum $
  else if ( sw_lun[20] eq 1 ) then $
     readu, LUN[20], dum

   XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
   if sw_debug eq 1 then begin
      print, 'H+ XIONN_m3=',jth,XIONN_m3[jth,60,0]

;dbg20151116 detect nan
      indicesN = where( dum eq !values.f_nan, countN )
      indicesZ = where( dum eq 0.0, countZ )
      if countN ne 0 then begin
         print, 'h+ Nan',countN,indicesN
         STOP
      endif
      if countZ ne 0 then begin
         print, 'h+ Zero',countZ,indicesZ
         STOP
      endif
   endif                        ;sw_debug eq 1 then  begin
endif                           ;( sw_lun

if ( sw_dif eq 0 ) then begin
   if ( sw_lun[8] eq 1 ) or (sw_lun[20] eq 1) then begin
      jth=2                     ;he+

      if ( sw_lun[8] eq 1 ) then $
         readu, LUN[8], dum $
      else if ( sw_lun[20] eq 1 ) then $
         readu, LUN[20], dum

      XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
      if sw_debug eq 1 then  begin
         print, 'He+ XIONN_m3=',jth,XIONN_m3[jth,60,0]

;dbg20151116 detect nan
         indicesN = where( dum eq !values.f_nan, countN )
         indicesZ = where( dum eq 0.0, countZ )
         if countN ne 0 then begin
            print, 'he+ Nan',countN,indicesN
            STOP
         endif
         if countZ ne 0 then begin
            print, 'he+ Zero',countZ,indicesZ
            STOP
         endif
      endif                     ;sw_debug eq 1 then  begin
endif                           ;   if ( sw_lun[8] eq 1 ) then begin



   if ( sw_lun[9] eq 1 ) or (sw_lun[20] eq 1) then begin
     jth=3                        ;n+


     if ( sw_lun[9] eq 1 ) then $
        readu, LUN[9], dum $
     else if ( sw_lun[20] eq 1 ) then $
        readu, LUN[20], dum

     XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
     if sw_debug eq 1 then  begin
        print, 'n+ XIONN_m3=',jth,XIONN_m3[jth,60,0]

;dbg20151116 detect nan
        indicesN = where( dum eq !values.f_nan, countN )
        indicesZ = where( dum eq 0.0, countZ )
        if countN ne 0 then begin
           print, 'n+ Nan',countN,indicesN
           STOP
        endif
        if countZ ne 0 then begin
           print, 'n+ Zero',countZ,indicesZ
           STOP
        endif
     endif                      ;sw_debug eq 1 then  begin
  endif                         ;( sw_lun

   if ( sw_lun[10] eq 1 ) or (sw_lun[20] eq 1) then begin
     jth=4                        ;no+

     if ( sw_lun[10] eq 1 ) then $
        readu, LUN[10], dum $
     else if ( sw_lun[20] eq 1 ) then $
        readu, LUN[20], dum

;dbg print, 'debug dum =',dum[60,0]
     XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
     if sw_debug eq 1 then  begin
        print, 'no+ XIONN_m3=',jth,XIONN_m3[jth,60,0]

;dbg20151116 detect nan
        indicesN = where( dum eq !values.f_nan, countN )
        indicesZ = where( dum eq 0.0, countZ )
        if countN ne 0 then begin
           print, '!STOP! INVALID NO+ Nan',countN,indicesN
;dbg20161128 stop commented out
           STOP
        endif
        if countZ ne 0 then begin
           print, '!STOP! INVALID NO+ Zero',countZ,indicesZ
;dbg20161128 stop commented out
           STOP
        endif
     endif                      ;sw_debug eq 1 then  begin
   endif ;( sw_lun

   if ( sw_lun[11] eq 1 ) or (sw_lun[20] eq 1) then begin
     jth=5                        ;o2+

     if ( sw_lun[11] eq 1 ) then $
        readu, LUN[11], dum $
     else if ( sw_lun[20] eq 1 ) then $
        readu, LUN[20], dum

     XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
     if sw_debug eq 1 then  begin
        print, 'O2+ XIONN_m3=',jth,XIONN_m3[jth,60,0]

;dbg20151116 detect nan
        indicesN = where( dum eq !values.f_nan, countN )
        indicesZ = where( dum eq 0.0, countZ )
        if countN ne 0 then begin
           print, 'o2+ Nan',countN,indicesN
           STOP
        endif
        if countZ ne 0 then begin
           print, 'o2+ Zero',countZ,indicesZ
           STOP
        endif
     endif                      ;sw_debug eq 1 then  begin
  endif                         ;( sw_lun

   if ( sw_lun[12] eq 1 )  or (sw_lun[20] eq 1) then begin
     jth=6                        ;n2+

     if ( sw_lun[12] eq 1 ) then $
        readu, LUN[12], dum $
     else if ( sw_lun[20] eq 1 ) then $
        readu, LUN[20], dum

     XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
     if sw_debug eq 1 then  begin
        print, 'n2+ XIONN_m3=',jth,XIONN_m3[jth,60,0]

;dbg20151116 detect nan
        indicesN = where( dum eq !values.f_nan, countN )
        indicesZ = where( dum eq 0.0, countZ )
        if countN ne 0 then begin
           print, 'n2+ Nan',countN,indicesN
           STOP
        endif
        if countZ ne 0 then begin
           print, 'n2+ Zero',countZ,indicesZ
           STOP
        endif
     endif                      ;sw_debug eq 1 then  begin
  endif                         ;( sw_lun

   if ( sw_lun[13] eq 1 ) or (sw_lun[20] eq 1) then begin
     jth=7                        ;o+(2D)

     if ( sw_lun[13] eq 1 ) then $
        readu, LUN[13], dum $
     else if ( sw_lun[20] eq 1 ) then $
        readu, LUN[20], dum

     XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
     if sw_debug eq 1 then  begin
        print, 'o+(2D) XIONN_m3=',jth,XIONN_m3[jth,60,0]

;dbg20151116 detect nan
        indicesN = where( dum eq !values.f_nan, countN )
        indicesZ = where( dum eq 0.0, countZ )
        if countN ne 0 then begin
           print, 'o+(2D) Nan',countN,indicesN
           STOP
        endif
        if countZ ne 0 then begin
           print, 'o+(2D) Zero',countZ,indicesZ
           STOP
        endif
     endif                      ;sw_debug eq 1 then  begin
  endif                         ;( sw_lun

   if ( sw_lun[14] eq 1 )  or (sw_lun[20] eq 1)  then begin
     jth=8                        ;o+(2P)

     if ( sw_lun[14] eq 1 )  then $
        readu, LUN[14], dum $
     else if ( sw_lun[20] eq 1 ) then $
        readu, LUN[20], dum

     XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
     if sw_debug eq 1 then  begin
        print, 'o+(2P) XIONN_m3=',jth,XIONN_m3[jth,60,0]

;dbg20151116 detect nan
        indicesN = where( dum eq !values.f_nan, countN )
        indicesZ = where( dum eq 0.0, countZ )
        if countN ne 0 then begin
           print, 'o+(2P) Nan',countN,indicesN
           STOP
        endif
        if countZ ne 0 then begin
           print, 'o+(2P) Zero',countZ,indicesZ
           STOP
        endif
     endif                      ;sw_debug eq 1 then  begin
  endif                         ;( sw_lun
endif ;( sw_dif eq 0 ) then begin

;Te:
if ( sw_lun[3] eq 1 ) then begin
   jth=3-1
   readu, LUN[3], dum
   TE_TI_k[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
   if sw_debug eq 1 then  print, 'TE_TI_k=',TE_TI_k[jth,0,0]
endif


;Ti:
if ( sw_lun[7] eq 1 ) then begin
   jth=1-1
   readu, LUN[7], dum
   TE_TI_k[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
   if sw_debug eq 1 then  print, 'TE_TI_k=',TE_TI_k[jth,0,0]
endif

if ( sw_dif eq 0 ) then begin
   if ( sw_lun[4] eq 1 ) then begin
      jth=0                     ;Vo+
      readu, LUN[4], dum
      XIONV_ms1[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
      if sw_debug eq 1 then  print, 'XIONV_ms1=',jth,XIONV_ms1[jth,60,0]
   endif ;( sw_lun[9] eq 1 ) then begin
endif

;sza
;if sw_debug eq 1 then  print,'sza sw_lun18', sw_lun[18]
;print, 'size',size(dum)
   if ( sw_lun[18] eq 1 ) then begin
      readu, LUN[18], dum
if sw_debug eq 1 then  print, 'dum sza_rad MAX=',MAX(dum)*180./!PI, MIN(dum)*180./!PI
      sza_rad[0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
      ;if sw_debug eq 1 then  $
;print, 'sza_rad=',sza_rad[60,0]*180./!PI
   endif ;( sw_lun[18] eq 1 ) then begin



if ( sw_dif eq 0 ) then begin
   dum=fltarr(NLP,NMP)
   if ( sw_lun[5] eq 1 ) then begin
      readu, LUN[5], dum
      jth=2
      for mp=0,nmp-1 do begin
;VEXBup(mp,lp)
         if  sw_debug eq 1 then  print,'VEXBup  mp=', mp
         VEXB[mp,0:nlp-1,jth]=dum[0:nlp-1,mp] ;+upward
      endfor                               ;mp
      if sw_debug eq 1 then   $
         print, 'VEXBup_ms1=',$
;VEXB[0,130,jth],VEXB[0,120,jth] ;low
VEXB[*,72,jth]
   endif                        ;( sw_lun[9] eq 1 ) then begin
;
   if ( sw_lun[15] eq 1 ) then begin
;VEXBe(mp,lp)
      readu, LUN[15], dum
jth=0
for mp=0,nmp-1 do begin
VEXB[mp,0:nlp-1,jth]=dum[0:nlp-1,mp] ;+eastward
endfor
      if sw_debug eq 1 then   $
;         print, 'VEXBe_ms1=',VEXB[0,130,jth],VEXB[0,120,jth]
;         ;original new grid 201207
         print, 'VEXBe_ms1=',VEXB[0,65,jth],VEXB[0,55,jth] ;2xdyn grid
   endif                        ;( sw_lun[9] eq 1 ) then begin
;
   if ( sw_lun[16] eq 1 ) then begin
;VEXBth(mp,lp)
      readu, LUN[16], dum
jth=1
      for mp=0,nmp-1 do begin
         VEXB[mp,0:nlp-1,jth]=dum[0:nlp-1,mp] ;+southward
      endfor 
      if sw_debug eq 1 then   $
         print, 'VEXBth_ms1=' $
;,VEXB[0,130,jth],VEXB[0,120,jth] ;original new grid 201207
,VEXB[0,65,jth],VEXB[0,55,jth] ;2xdyn grid
   endif                        ;( sw_lun[9] eq 1 ) then begin
endif

;sunlons
if ( sw_lun[17] eq 1 ) then begin
  readf, LUN[17],sunlons1
  print,'sunlons1=', sunlons1
endif
  

if ( sw_3DJ eq 1 ) then begin
     readf, LUN[6], je_3d
;if sw_debug eq 1 then  
print, 'je_3d[unit???]=',je_3d[0,0:20,0]
endif

if ( sw_hr eq 1 ) then begin
mptmp=0L
hrtmp=fltarr(NPTS2D)

for jth=1,7 do begin
;hr-jth
     readf, LUN[jth+5], mptmp
     readf, LUN[jth+5], hrtmp ;(1) !!PGR elec_ion
;jth=1-1
hrate[jth-1,1-1:NPTS2D-1,mptmp-1]=hrtmp[1-1:NPTS2D-1]
endfor

endif

;;;;dbg transport
;jth=1L
;mp=0L
;dum=fltarr(NPTS2D)
;     readf, LUN[4], dum
;XIONN_m3[jth,0:NPTS2D-1,mp]=dum[0:NPTS2D-1]
;if sw_debug eq 1 then  print, 'XIONN_m3=',jth,XIONN_m3[jth,60,mp]



END ;PRO read_bin
