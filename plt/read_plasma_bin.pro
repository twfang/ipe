pro read_plasma_bin,LUN,UT_hr, XIONN_m3,XIONV_ms1,TE_TI_k,VEXB,sw_debug $
,sw_3DJ,je_3d,sw_hr,hrate, sw_dif

UT_sec=0L
record_number=0L
; read plasma1
     readf, LUN[0], record_number, UT_sec
 UT_hr = UT_sec /3600.
print,' rec#',record_number,' UThr',UT_hr, UT_sec

size_result=size(XIONN_m3)
if ( sw_debug eq 1 ) then  print,size_result
NPTS2D=size_result[2]
NMP=size_result[3]
dum=fltarr(NPTS2D,NMP)
jth=0 ;o+
     readu, LUN[2], dum
;dbg print ,'check o+',MAX(dum),MIN(dum)
XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
if sw_debug eq 1 then  print, 'o+ XIONN_m3=',jth,XIONN_m3[jth,60,0]

jth=1 ;h+
     readu, LUN[6], dum
XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
if sw_debug eq 1 then  print, 'H+ XIONN_m3=',jth,XIONN_m3[jth,60,0]

if ( sw_dif eq 0 ) then begin
jth=2 ;he+
     readu, LUN[8], dum
XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
if sw_debug eq 1 then  print, 'He+ XIONN_m3=',jth,XIONN_m3[jth,60,0]
jth=3 ;n+
     readu, LUN[9], dum
XIONN_m3[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
if sw_debug eq 1 then  print, 'N+ XIONN_m3=',jth,XIONN_m3[jth,60,0]
endif

;Te:
jth=3-1
     readu, LUN[3], dum
TE_TI_k[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
if sw_debug eq 1 then  print, 'TE_TI_k=',TE_TI_k[jth,0,0]

;Ti:
jth=1-1
     readu, LUN[7], dum
TE_TI_k[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
if sw_debug eq 1 then  print, 'TE_TI_k=',TE_TI_k[jth,0,0]

if ( sw_dif eq 0 ) then begin
jth=0;Vo+
     readu, LUN[4], dum
XIONV_ms1[jth,0:NPTS2D-1,0:NMP-1]=dum[0:NPTS2D-1,0:NMP-1]
if sw_debug eq 1 then  print, 'XIONV_ms1=',jth,XIONV_ms1[jth,60,0]
endif

if ( sw_dif eq 0 ) then begin
;VEXB(mp,lp)
     readu, LUN[5], VEXB
;if sw_debug eq 1 then  print, 'VEXB_ms1=',VEXB[0,100,0],VEXB[1,100,0]
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
