;date: 20150616
;version: only ed1 for CJ
pro rd_ed12_20150616


; choose jicamarca
   mp0=0L
   lp0=130-1L                    ;-9.05[deg] mp=0;jicamarca: ht=254.107km


; read fort.7000 for Be and apexE
input_DIR0='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r311.1tmp/trunk/run/ipe_S_1539'
openr, LUN7, input_DIR0+'/fort.7000', /GET_LUN
readf, lun7, mp
print, mp
readf, lun7, Be3
print, Be3
;readf, lun7, apexE1
;print, apexE1
;readf, lun7, apexE2
;print, apexE2

;close file
free_lun,lun7

input_DIR='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/reu/tmp20130703reu/trunk/run/ipe_640_18702'

openr, LUN6, input_DIR+'/ut_rec', /GET_LUN ;utime
openr, LUN8, input_DIR+'/fort.2008', /GET_LUN ;ed1_90
openr, LUN9, input_DIR+'/fort.2009', /GET_LUN ;ed2_90


formatE='(20E12.4)'
nmp=80L
nlp=170L
ed190=fltarr(nmp,nlp*2) ;[mV/m]
ed290=fltarr(nmp,nlp*2) ;[mV/m]


nmin=1-1L
nmax=97-1L;97L
vExBUp=fltarr(nmax+1)
vExBe =fltarr(nmax+1)
lthr=fltarr(nmax+1)
for n_read=nmin,nmax do begin

;read ed190
   edum=fltarr(2,nlp,nmp)
   readf, LUN8, edum,  FORMAT=formatE ;ed190
   for mp=0,nmp-1 do begin
      ii=-1
      for lp=0,nlp-1 do begin
         ii=ii+1
         ed190[mp,ii] = edum[1,lp,mp] ;SH
      endfor                          ;lp
      for lp=nlp-1,0,-1 do begin
         ii=ii+1
         ed190[mp,ii] = edum[0,lp,mp] ;NH
      endfor                          ;lp
   endfor                             ;mp

;read ed290
   edum=fltarr(2,nlp,nmp)
   readf, LUN9, edum,  FORMAT=formatE ;ed290
   for mp=0,nmp-1 do begin
      ii=-1
      for lp=0,nlp-1 do begin
         ii=ii+1
         ed290[mp,ii] = edum[1,lp,mp]*(-1.) ;SH
      endfor                          ;lp
      for lp=nlp-1,0,-1 do begin
         ii=ii+1
         ed290[mp,ii] = edum[0,lp,mp]*(-1.) ;NH
      endfor                          ;lp
   endfor                             ;mp
 


   v_e1 =   ed290[mp0,lp0]*1.E-3 / Be3 ;!(4.18) +mag-east 
   v_e2 = - ed190[mp0,lp0]*1.E-3 / Be3 ;!(4.19) +down/equatorward
   
   vExBUp[n_read] = v_e2 * (-1.)  ;down-->upward
   vExBe[ n_read] = v_e1          ;
   print, 'n_read=', n_read,' ed1[mV/m]=',ed190[mp0,lp0],' ed2[mV/m]=',ed290[mp0,lp0]
   print,' vEXBup[m/s]=', vExBup[n_read],' vEXBe[m/s]=', vExBe[n_read]


;read Universal Time (UT)
rec=0L
utime=0L
ut0=168.
   readf, LUN6, rec, utime           ;in seconds
   uthr = utime / 3600. - ut0
;calculate Local Time (LT)
   lthr[n_read] = uthr - 5.
   print, 'utime=',utime,' UT hr=', uthr,' LT=',lthr[n_read]
   
endfor                          ;n_read=nmin,nmax

;close file
free_lun,lun8


;plot EXB drift as a function of LT
;set up a display window
DEVICE, RETAIN=2, DECOMPOSED=0
WINDOW,1,XSIZE=1000,YSIZE=1000
loadct, 0  ;black and white

;plot vEXB upward component
plot,lthr,vExBup $
, yrange=[-90. , +90. ],  ystyle=1  $
, xrange=[ MIN(lthr), MIN(lthr) ],  xstyle=1  $
,XTITLE = 'LT[hr] ', YTITLE = 'VEXB[m/s]' $
, charsize =2.5, charthick=2.5 $
, TITLE=' ' $
, linestyle=0 $
, thick=2.0 $
, color=250.

;overplot the eastward drift 
oplot,lthr,vExBe $
, linestyle=4 $
, thick=2.0 $
, color=250.

;save figure in a file
output_png, 'vexb.png'


; output the values to ascii files:
openw,lun5,'output_lthr', /GET_LUN
printf, lun5,lthr
free_lun,lun5

openw,lun5,'output_exbup', /GET_LUN
printf, lun5,vExBup
free_lun,lun5

openw,lun5,'output_exbe', /GET_LUN
printf, lun5,vExBe
free_lun,lun5

print, 'program rd_ed12 finished!'
end ;pro rd_ed12
