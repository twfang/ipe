pro get_tec_time_variation ;
;,time_string $
;,tec4,nmf24,hmf24 $
;,tec5,nmf25,hmf25 $
;,tec6,nmf26,hmf26




TEST1=$
;'1452876517' ;jid20160112 v7 16.5ut
'1452876646' ;jid20160112 v8 16ut
;'1452823230' ori
;time_string='137.0000UT'

path='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/fig/tmp20150730/'+TEST1+'_ipe_theia_intel_parallel_80/'

uthr0=136.0;120.00

ntmin=0
ntmax=3+2;68;97
utimesv=fltarr(ntmax)
tecsv=fltarr(ntmax)
nmf2sv=fltarr(ntmax)
hmf2sv=fltarr(ntmax)
time_step=15./60.

for ntime=ntmin,ntmax-1 do begin
   uthr = uthr0 + (ntime)*time_step
   print, ntime, uthr
   if ( uthr lt 1.0E+01 ) then $
      time_string='00'+STRTRIM(STRING( uthr, FORMAT='(F8.4)'),1)+'UT'  $
   else if ( uthr lt 1.0E+02 ) then $
      time_string='0'+STRTRIM(STRING( uthr, FORMAT='(F8.4)'),1)+'UT'  $
   else $                       ;if ( uthr ge 1.0E+02 ) then $
      time_string=STRTRIM(STRING( uthr, FORMAT='(F8.4)'),1)+'UT'
   
   print, time_string
   restore ,path+'TEC_tmp20150730_ipe_80_'+TEST1+'_'+time_string+'.sav'
   
;get i for mlat=-30
   iplt=30
   utimesv[ntime]=uthr - 120.00
   tecsv[ntime]=tec0[iplt]
   nmf2sv[ntime]=nmf20[iplt]
   hmf2sv[ntime]=hmf2[iplt]

endfor

 !P.BACKGROUND=255

plot, utimesv,nmf2sv $
, TITLE='nmf2: '+TEST1 $ ;+FileID $
, XTITLE = 'UT [hours]',    YTITLE = 'NmF2'  $
, yrange=[ 0., 2.7 ], ystyle=1  $ ;20141119
, xrange=[ 0., +17. ], xstyle=1  $
,color=text_color, THICK =line_thick, charsize=char_size, charthick=char_thick

output_png, path+'nmf2TimeVariation.png'
save, /VARIABLES, filename=path+'TEC_TimeVariation_tmp20150730_ipe_80_'+TEST1+'.sav'

print, 'get_TecTimeVariation finished '
end


