;date: 20151009
; purpose: read f107 and plot them
; for GRL reviewer
pro rd_f107 
   
sw_output2file=1
if ( sw_output2file eq 1 ) then flnm_png='f107Kp.png'

  nday0=341L
nmax=31L
f107=fltarr(31)
nday=intarr(31)

;---(1)read 2012
input_file='f107adj2012.txt'
openr, lun, input_file, /get_lun

str_tmp='==============================================================================='
for n = 0, 4 do begin
   readf, lun, str_tmp, FORMAT='(A82)'
print,str_tmp
endfor ;n


iday=0L
idum1='123456'
idum2=0L
idum3=0L
idum4=0L
idum5=0L
idum6=0L
idum7=0L
idum8=0L
idum9=0L
idum10=0L
idum11=0L
idum12=0L

ii=-1L
for n = 1, 31 do begin

   readf, lun, iday,$
          idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11,idum12,$
          FORMAT='(i5,A6,11i6)'

   print,iday,$
         idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11,idum12,$
         FORMAT='(i5,A6,11i6)'

   if ( (n mod 5) eq 0 ) then begin
      readf, lun, str_tmp, FORMAT='(A82)'
      print,n,'blank',str_tmp
   endif
   
   if n lt 6 then CONTINUE      ;go to next 

   ii = ii + 1
   ;take Dec data
   f107[ii] = FIX(idum12) * 0.1
   nday[ii] = ii + nday0
                                ;d print,n,' ii=',ii, ' f107=',(idum12 * 0.1)


endfor; n = 0, nmax do begin
free_lun, lun
;d print,ii,f107
;d plot, f107


;---(2)read 2013 f107
input_file='f107adj2013.txt'
openr, lun, input_file, /get_lun

str_tmp='==============================================================================='
for n = 0, 4 do begin
   readf, lun, str_tmp, FORMAT='(A82)'
print,str_tmp
endfor ;n



idum1=0L
for n = 1, 5 do begin

   readf, lun, iday,$
          idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11,idum12,$
          FORMAT='(i5,12i6)'

   print,iday,$
         idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11,idum12,$
         FORMAT='(i5,12i6)'

;   if ( (n mod 5) eq 0 ) then begin
;      readf, lun, str_tmp, FORMAT='(A82)'
;      print,n,'blank',str_tmp
;   endif
   
   ii = ii + 1
   ;take Jan data
   f107[ii]= FIX(idum1) * 0.1
   nday[ii] = ii + nday0
   print,n,' ii=',ii, ' f107=',(idum1 * 0.1)


endfor; n = 0, nmax do begin
free_lun, lun

print,ii,f107

;---read kp and ap
  nmax1=2L;tmp 37L
input_file='Kpap201212_201301.txt'
openr, lun, input_file, /get_lun

str_tmp='YYYYMMDD Kp[8]           Sum ap[8]                   Ap'
readf, lun, str_tmp, FORMAT='(A55)'
print, str_tmp


;ap=LONarr(32)
kp=LONarr(32*8)

str_tmp='20121201 1 1 0+0+1-1-2+2+ 9-  4  4  2  2  3  3  9  9'
idum0='20121201'
idum1=0L
ii=-1L
for n = 0, nmax1-1 do begin

readf, lun, idum0 $
,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8 $
, FORMAT='(i8,8i2)'

print ,idum0,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8 , FORMAT='(i9,8i3)'

if n lt 5 then CONTINUE ;go to next 


ii = ii + 1
kp[ii]= idum1
ii = ii + 1
kp[ii]= idum2
ii = ii + 1
kp[ii]= idum3
ii = ii + 1
kp[ii]= idum4
ii = ii + 1
kp[ii]= idum5
ii = ii + 1
kp[ii]= idum6
ii = ii + 1
kp[ii]= idum7
ii = ii + 1
kp[ii]= idum8







print,n,' ii=',ii, ' kp=',idum[0], idum[1]



endfor; n = 0, nmax1 do begin
STOP ;tmp
free_lun,lun

print,MIN(ap),MAX(ap)
;d plot, ap

print, 'pro rd_ap finished'


;---plotting     
iwindow=1
size_window=1000
device, retain=2, decomposed=0
window,iwindow, xsize=size_window,ysize=size_window
!P.BACKGROUND=255



char_size=3.0
char_thick=3.0
line_thick=4.0

x_min=341;min(nday)
x_max=371;max(nday)
y_min=90.;min(f107)
y_max=140.;max(f107)
print, x_min,x_max
print, y_min,y_max

loadct,0
colory=0.
plot,nday,f107 $
     ,TITLE = 'F!D10.7!N & A!DP' $
     ,XTITLE = 'DAY OF YEAR' $
     ,YTITLE = 'F!D10.7' $
     ,xstyle=1,  ystyle=9, xmargin=[8,8], ymargin=[4,4] $
, xrange=[x_min, x_max ] $
, yrange=[y_min ,y_max ] $
, color=colory           $
     , charsize=char_size, charthick=char_thick $
     , thick = line_thick
;, position=[0.2,0.2,0.8,0.8] $
;, /NODATA
     

;     oplot,nday,f107 $
;           , color=colory


; Draw the top x-axis, supplying labels, etc.
; Make the characters smaller so they will fit:
;MONTH=['DEC','JAN']
;AXIS, XAXIS=1, XTICKS=1, XTICKV=NDAY, XTICKN=MONTH, $
;      XTITLE='Month', XCHARSIZE = 0.7


                                ;Draw the right y-axis

     AXIS,/SAVE, YAXIS=1, YRANGE=( (!Y.CRANGE-90.)*.25 ), YSTYLE = 1 $
          , YTITLE = 'A!DP' $
          , color=colory  $
          , charsize=char_size, charthick=char_thick
;          ;, YSTYLE=1, $

     oplot,nday,ap $
          , linestyle=5 $
           , color=colory $
                , thick = line_thick


if ( sw_output2file eq 1 ) then output_png, flnm_png

print, 'pro rd_f107 finished'
end ;pro rd_f107
