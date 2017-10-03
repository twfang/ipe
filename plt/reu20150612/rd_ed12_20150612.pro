;date: 20150612
;version: only ed1 for CJ
pro rd_ed12


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
readf, lun7, apexE1
print, apexE1
readf, lun7, apexE2
print, apexE2

;close file
free_lun,lun7

input_DIR='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/reu/tmp20130703reu/trunk/run/ipe_640_18702'

openr, LUN6, input_DIR+'/ut_rec', /GET_LUN ;utime
openr, LUN8, input_DIR+'/fort.2008', /GET_LUN ;ed1_90


formatE='(20E12.4)'
nmp=80L
nlp=170L
ed190=fltarr(nmp,nlp*2) ;[mV/m]


nmin=1-1L
nmax=97-1L;97L
vExBGeoUp=fltarr(nmax+1)
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

 


;!            v_e(1) =   Ed2_90(1,lp,mp) / Be3(lp,mp) !(4.18) +mag-east(d1?) 
   v_e2 = - ed190[mp0,lp0]*1.E-3 / Be3 ;!(4.19) +down/equatorward(d2?)
   
   vExBGeoUp[n_read] = $
                                ;             (v_e1*apexE1) + $       ;& 
      (v_e2*apexE2)
   print, 'n_read=', n_read,' ed1[mV/m]=',ed190[mp0,lp0],' vEXB[m/s]=', vExBGeoUp[n_read]


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


openw,lun5,'output_lthr', /GET_LUN
printf, lun5,lthr
free_lun,lun5

openw,lun5,'output_exb', /GET_LUN
printf, lun5,vExBGeoUp
free_lun,lun5

print, 'program rd_ed12 finished!'
end ;pro rd_ed12
