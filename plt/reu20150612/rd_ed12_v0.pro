;date: 20150612
;v0: both ed1 & ed2
pro rd_ed12


input_DIR='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/reu/tmp20130703reu/trunk/run/ipe_640_18702'
LUN8=0L
;LUN9=0L
openr, LUN8, input_DIR+'/fort.2008', /GET_LUN ;ed1_90
;openr, LUN9, input_DIR+'/fort.2009', /GET_LUN ;ed2_90

formatE='(20E12.4)'
nmp=80L
nlp=170L
ed190=fltarr(nmp,nlp*2)
;ed290=fltarr(nmp,nlp*2)

nmin=1L
nmax=97L
for n_read=nmin,nmax do begin
   print, 'n_read=', n_read
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
   
;  readf, LUN9, edum,  FORMAT=formatE   ;ed290
;  for mp=0,nmp-1 do begin
;    ii=-1
;    for lp=0,nlp-1 do begin
;      ii=ii+1
;      ed290[mp,ii] = edum[1,lp,mp];SH
;    endfor ;lp
;    for lp=nlp-1,0,-1 do begin
;      ii=ii+1
;      ed290[mp,ii] = edum[0,lp,mp];NH
;    endfor ;lp
;  endfor ;mp

   mp=0L
   lp=130-1L                    ;-9.05[deg] mp=0;jicamarca: ht=254.107km
   print,'ed1=',ed190[mp,lp]
;print,'ed2=',ed290[mp,lp]
   
endfor                          ;n_read=nmin,nmax

print, 'program rd_ed12 finished!'
end ;pro rd_ed12
