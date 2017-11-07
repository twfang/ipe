;20140117: quick plot version
;20131219: hinotori: Vartype in ctr_lon_lat needs to be zero to be able to output Ne. output files are saved in ~/wamns/hinotori/. sw_plot_contour=0 to run faster by saving plotting time.
pro ctr_lon_lat_quick $
,JMIN_IN,JMAX_IS,Z_km,mlat_deg $
, XIONN_m3, TE_TI_k $
, XIONV_ms1 $
,UT_hr, plot_DIR $
,n_read $
,sw_output2file $
,glon_deg,glat_deg,sw_frame,fac_window ,TEST $
,sw_debug $
;20131209: output to ascii file
,sw_output2file_ascii,luntmpN,luntmpS,ncount $
, Vn_ms1,tn_k,on_m3,n2n_m3,o2n_m3 $
, n_plt_max,input_DIR0 $
, ht_plot,rundir $
, VarType_max, VarType_min, VarType_step $
, n_plt_min ;=0L
print, 'ht_plot=', ht_plot
;print, 'VarType_max', VarType_max, 'VarType_min', VarType_min
;print,  'VarType_step',  VarType_step


print,'inside ctr_lon_lat_quick: check on_m3', MAX( on_m3 ), MIN( on_m3 )


;nm20161220 profile validation
sw_debug_prfl=0L

 sw_output_nmf2=0 ;NH:1; N&S:2 ; no nmf2 output:0
 sw_output_nmf2noon=0
;dbg20140806
;d print,'ctr_lon_lat_quicj:check vn_ms1', MAX( vn_ms1(2-1,*,*) ), MIN( vn_ms1(2-1,*,*) )

; X-axis range
;1:0<lon<360
;0:-180<lon<+180
xmax360=getenv('xmax360')
;note20131209: sw OFF only when sw_output2file_ascii=1 to run faster!!
sw_plot_contour=1
dlt=360./80./360.*24.
ltimemin=12. - dlt*.5
ltimemax=12. + dlt*.5
;print,' ltimemin', ltimemin, 'ltimemax', ltimemax
mlatmax=63.
factor=1.0E-11 ;E-12 ;for density 
; remember to modify zmin/max

;debug20140108
fac_wind=1.0;1.0E-2

;VarTitle=['Ne','Te','Ti','o+','NO+','O2+','Vpar_o+','NmF2','HmF2']
VarTitle=['Ne','Te','Ti','o+','TN','UN'$
;,'VN'$
,'On_m3'$  ;[O] [m-3]
,'NmF2','HmF2'] ;20140108
;VarTitle=['Ne','Te','Ti','N(NO+)','Vpar_o+']
;VarTitle=['Ne','Te','Ti','N(O2+)','Vpar_o+']
;VarTitle=['Ne','Te','Ti','N(N2+)','Vpar_o+']

unit=['[10^11 m-3]','K','K','[m-3]',$
;'[m-3]'$
'[K]'$
,'[m/s]'$
;,'[m/s]'$
,'[log10 m-3]'$
,'[10^11 m-3]','[km]'] ;[nA/m2]'];[nA/m2]' ;20140108

;20140203: remember Vartype loop cannot be used for quick plot version!!!
for VarType=VarType_min, VarType_max, VarType_step do begin

if ( n_read eq n_plt_min ) then  print,'VarType=',VarType




sw_range=1L
nano=1.0E-9


 sw_plot_grid=0 ;1:retangular, 2:polar
; get n_read_max, ny_max
size_result = SIZE(JMIN_IN)
if ( sw_debug eq 1 ) and ( n_read eq n_plt_min ) then  print,'NLP',size_result
NLP = size_result[1]
ny_max =NLP*2
if  ( sw_debug eq 1 ) and ( n_read eq n_plt_min ) then $
   print,'NLP',NLP,'ny_max',ny_max

size_result = SIZE(XIONN_m3) ;je_3d)
if  ( sw_debug eq 1 ) and ( n_read eq n_plt_min ) then  print,'NMP',size_result
ISPEC=size_result[1]
NPTS2D=size_result[2]
NMP = size_result[3]
nx_max =size_result[3]    ;NMP=80
if  ( sw_debug eq 1 ) and ( n_read eq n_plt_min ) then $
  print,'NMP=',NMP,' nx_max=',nx_max,' ISPEC=',ispec,' NPTS2D=',NPTS2D


mlon_deg_min=0. ;default
dmlon=4.50 ;deg
mlon_deg=findgen(NMP)*dmlon
for mp=0,nmp-1 do begin
  if ( mlon_deg[mp] lt mlon_deg_min ) then  mlon_deg[mp]=mlon_deg[mp]+360.
endfor

plot_zz=fltarr(nx_max,ny_max)
plot_yy=fltarr(nx_max,ny_max)
plot_xx=fltarr(nx_max,ny_max)
;nm20140926 wind output
plot_zz1=fltarr(nx_max,ny_max)



for mp=0,NMP-1 do begin
   for lp=0,NLP-1 do begin


  in = JMIN_IN[lp]-1L
  is = JMAX_IS[lp]-1L
  midpoint = IN + ( IS - IN )/2
;NH  
  nel = fltarr(NPTS2D)*0.0000
  istep=+1
  for i=in,midpoint, istep  do begin

   if ( VarType lt 7 ) then begin
;dbg print,in,'!dbg',i,' mp=',mp,' lp=',lp,z_km[i],ht_plot,z_km[i+istep],on_m3[i,mp]


if ( z_km[i] le ht_plot ) AND (z_km[i+istep] gt ht_plot ) then begin


;lp=68
;dbg if mp eq 64 AND lp ge 58 AND lp le 78 then $
;dbg print,'NH!dbg',i,' mp=',mp,' lp=',lp,XIONN_m3[0:8,i,mp],z_km[i]




      if ( VarType eq 0 ) then $ 
        for jth=0,ISPEC-1 do  plot_zz[mp,lp]=plot_zz[mp,lp]+XIONN_m3[jth,i,mp] * factor $
      else if ( VarType eq 1 ) then $
        plot_zz[mp,lp] = TE_TI_k[3-1,i,mp] $;[K] ;Te
      else if ( VarType eq 2 ) then $
        plot_zz[mp,lp] = TE_TI_k[1-1,i,mp] $;[K] ;Ti
      else if ( VarType eq 3 ) then $
        plot_zz[mp,lp] = XIONN_m3[1-1,i,mp]*factor  $ ;o+
      else if ( VarType eq 4 ) then $
        plot_zz[mp,lp] = $
;XIONN_m3[4,i,mp]*factor  $ ;no+
tn_k[i,mp] $;[k] ;tn
      else if ( VarType eq 5 ) then $
        plot_zz[mp,lp] =$
Vn_ms1[2-1,i,mp]*fac_wind $;[m/s] ;positive northward ;20140108
; XIONN_m3[5,i,mp]*factor  $ ;o2+
      else if ( VarType eq 6 ) then begin ;$
        plot_zz[mp,lp] =$
alog10(on_m3[i,mp]) ;[m-3] ; oxygen density
;Vn_ms1[2-1,i,mp]*fac_wind ;[m/s] ;
; XIONV_ms1[1-1,i,mp] ;V//o+[m/s]

;if sw_debug eq 1 then 
print, i,mp, z_km[i],mlat_deg[i],on_m3[i,mp], plot_zz[mp,lp]
endif

 
      if ( sw_frame eq 0 ) then begin ;magnetic
         plot_yy[mp,lp] = mlat_deg[i]


;if (mp eq 0) and (lp eq 46) then begin
;print,'check mlat', mlat_deg[i]
;STOP
;endif

         plot_xx[mp,lp] = mlon_deg[mp]
         if (xmax360 eq 0) AND (mlon_deg[mp] ge +180.) then plot_xx[mp,lp]=mlon_deg[mp]-360.
      endif else if ( sw_frame eq 1 ) then begin ;geographic
         plot_yy[mp,lp] = glat_deg[i,mp]
         plot_xx[mp,lp] = glon_deg[i,mp]
         if (xmax360 eq 0 ) AND (glon_deg[i,mp] ge +180.) then plot_xx[mp,lp]=glon_deg[i,mp]-360.

;      endif else if ( sw_frame eq 2 ) then begin ;geographic
;         plot_yy[mp,lp] = mlat_deg[i]
;;
;
;;dbg20140923: 129 is jicamarca
;  ltime = ut_hr + glon_deg[i,mp]/15.
;  if ( ltime ge 24. ) then ltime = ( ltime MOD 24. )
;;extract near noon data
;if ( lp eq 129 ) then  begin
;  if ltime ge ltimemin and ltime le ltimemax then  $
;    print, 'NH: mp',mp,'lp',lp,'LT',ltime, ut_hr, mlat_deg[i],glon_deg[i,mp]
;  mp_output_noonprfl = mp
;  glon_output = glon_deg[i,mp]
;endif
;
;         plot_xx[mp,lp] = ltime

      endif ;( sw_frame eq 0 ) then begin ;magnetic



;20131209output to file:(1)NH
;20140808ipe 1st paper
       if ( sw_output2file_ascii eq 1 ) AND ( ABS( mlat_deg[i]) le mlatmax ) then begin
          ltime=UT_hr + glon_deg[i,mp]/15.
          if ( ltime ge 24. ) then ltime = ( ltime MOD 24. )
          if ( ltime ge ltimemin) AND ( ltime le ltimemax ) then begin
             ncount=ncount+1
             printf, luntmp, UT_hr, plot_zz[mp,lp], TE_TI_k[3-1,i,mp], TE_TI_k[1-1,i,mp],mlat_deg[i], glon_deg[i,mp], ltime 
          endif                 ;( ltime ge 11.) AND ( ltime le 15. ) then begin
       endif                    ;( sw_output2file_ascii eq 1 ) then 


      BREAK ;exit from the i loop 
    endif                       ;( z_km[i] le ht_plot ) AND ( z_km[i+istep] gt ht_plot ) then begin

;nm20130523: added nmf2/hmf2
  endif else begin  ;( VarType ge 7 ) then begin
     for jth=0,ISPEC-1 do  nel[i]=nel[i]+XIONN_m3[jth,i,mp]

     if ( i eq midpoint ) then begin
       result = MAX( nel,Max_Subscript ) ;find NmF2
;       print,'NH: NmF2', result,Max_Subscript, nel[Max_Subscript],Z_km[Max_Subscript],mp,lp
       if ( VarType eq 7 ) then $
         plot_zz[mp,lp] = nel[Max_Subscript] * factor $ ;NmF2
       else if ( VarType eq 8 ) then $
         plot_zz[mp,lp] = Z_km[Max_Subscript]    ;hmF2

       if ( sw_frame eq 0 ) then begin ;magnetic
          plot_yy[mp,lp] = mlat_deg[Max_Subscript]
          plot_xx[mp,lp] = mlon_deg[mp]
          if (xmax360 eq 0) AND (mlon_deg[mp] ge +180.) then plot_xx[mp,lp]=mlon_deg[mp]-360.

       endif else if ( sw_frame eq 1 ) then begin ;geographic
          plot_yy[mp,lp] = glat_deg[Max_Subscript,mp]
          plot_xx[mp,lp] = glon_deg[Max_Subscript,mp]
          if (xmax360 eq 0) AND (glon_deg[Max_Subscript,mp] ge +180.) then plot_xx[mp,lp]=glon_deg[Max_Subscript,mp]-360.

;dbg20140923 NH
       endif else if ( sw_frame eq 2 ) then begin ;noon in the middle
;dbg20140923: 129 is jicamarca
  ltime = ut_hr + glon_deg[Max_Subscript,mp]/15.
  if ( ltime ge 24. ) then ltime = ( ltime MOD 24. )
;extract near noon data
;if ( lp eq 129 ) then  begin
;   if ltime ge ltimemin and ltime le ltimemax then begin
;      print, 'NH: mp',mp,'lp',lp,'LT',ltime, ut_hr, mlat_deg[Max_Subscript],glon_deg[Max_Subscript,mp]
;      mp_output_noonprfl = mp
;      glon_output = glon_deg[Max_Subscript,mp]
;   endif
;endif

         plot_yy[mp,lp] = mlat_deg[Max_Subscript]
         plot_xx[mp,lp] = ltime
;nm20150225         plot_zz1[mp,lp]= Vn_ms1[2-1,Max_Subscript,mp]*fac_wind ;field aligned wind positive NORTHward



       endif ;( sw_frame eq 0 ) then begin ;magnetic
  

;dbg20140911validate wind NH
;dbg20140825 validate wind direction nmf2
;if ( mp eq 51 ) then print, mp,lp, lp
if ( sw_output_nmf2 ge 1 AND mp eq 51 AND lp eq 36 ) then begin
if ( n_read eq n_plt_min ) then begin
;   luntmpN=100
   flnmtmpN='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r319/trunk/run/'+rundir+'/nmf2NH.dat'
   openw,luntmpN,flnmtmpN, /GET_LUN
   print, 'NH nmf2 file created:',flnmtmpN
endif
;dbg20140825
     ltime=UT_hr + glon_deg[i,mp]/15.
     if ( ltime ge 24. ) then ltime = ( ltime MOD 24. )
print, mp,lp,UT_hr, plot_zz[mp,lp], Z_km[Max_Subscript], glat_deg[Max_Subscript,mp], glon_deg[Max_Subscript,mp], ltime
Printf, luntmpN, UT_hr, plot_zz[mp,lp], Z_km[Max_Subscript], glat_deg[Max_Subscript,mp], glon_deg[Max_Subscript,mp], ltime
endif ;( lp eq 47 AND mp eq 49 ) then begin


    endif;     if ( i eq midpoint ) then 
  endelse ;( VarType lt 7 ) then begin



  endfor ;i=in,midpoint, istep  do begin

;SH  
  lps = NLP-1 + (NLP-1-lp)
  nel = fltarr(NPTS2D)*0.0000
  istep=-1 
  for i=is,midpoint, istep  do begin
   
    if ( VarType lt 7 ) then begin



if ( z_km[i] le ht_plot ) AND ( z_km[i+istep] gt ht_plot ) then begin


;lps=270
;dbg if mp eq 64  AND lps ge 260 AND lps le 280 then begin
;dbg    print,'SH!dbg',i,' mp=',mp,' lp=',lp,' lps=',lps,XIONN_m3[0:8,i,mp],z_km[i]
;if mp eq 48 AND lps eq 270 then print,'!dbg',i,mp,lps,XIONN_m3[0:8,i,mp]
   ;STOP
;dbg endif
      if ( VarType eq 0 ) then $ 
        for jth=0,ISPEC-1 do  plot_zz[mp,lps]=plot_zz[mp,lps]+XIONN_m3[jth,i,mp]*factor $
      else if ( VarType eq 1 ) then $
        plot_zz[mp,lps] = TE_TI_k[3-1,i,mp] $;[K]
      else if ( VarType eq 2 ) then $
        plot_zz[mp,lps] = TE_TI_k[1-1,i,mp] $;[K]
      else if ( VarType eq 3 ) then $
        plot_zz[mp,lps] = XIONN_m3[1-1,i,mp]*factor $
      else if ( VarType eq 4 ) then $
        plot_zz[mp,lps] = $
;XIONN_m3[4,i,mp]*factor $
tn_k[i,mp] $;[k] ;Tn
      else if ( VarType eq 5 ) then $
        plot_zz[mp,lps] =$
Vn_ms1[1-1,i,mp]*fac_wind $;[m/s] ;eastward  ;20140108
; XIONN_m3[5,i,mp]*factor $
      else if ( VarType eq 6 ) then $
        plot_zz[mp,lps] = $
alog10(on_m3[i,mp]) ;[m-3] ;oxygen density
;Vn_ms1[2-1,i,mp]*fac_wind $;[m/s] ;northward  ;20140108
;XIONV_ms1[1-1,i,mp] ;V//o+[m/s]

if sw_debug_prfl eq 1 then begin
   mpPlot=68L
lpPlot=35L ;SH -49.79  218.11 NH  40.19  243.17
;lpPlot=0L; SH -74.69  133.07 NH  80.19  269.80
; lpPlot=120L                  ; SH -15.12  232.16 NH   6.87  236.12
   YrangeMax=1000.
;   YrangeMax=1060.
   if lpPlot eq 35L then $
      ;XrangeMin=1.E+4 $;1.E-10 $
      XrangeMin=1.E-10 $
   else if lpPlot eq 0L then $
      XrangeMin=1.E+4 $ ;1.E-15 $
   else if lpPlot eq 120L then begin
      XrangeMin=1.E+5
      YrangeMax=350.
   endif
  
;whichModel='WAM'
 whichModel='MSIS'
   if whichModel eq 'MSIS' then $
      n_plot=0L $               ;MSIS
   else if whichModel eq 'WAM' then $
      n_plot=1L                 ;WAM
   
   if n_read eq n_plot AND mp eq mpPlot AND lp eq lpPlot then begin

      ltmp=(ut_hr-120.+glon_deg[in,mp]/15.)
      print,mp,lp,' LT=',ltmp
      
      print, in,' mp=',mp,' lp=',lp,lps,mlat_deg[in],mlon_deg[mp],' glat=',glat_deg[in,mp],' glon=',glon_deg[in,mp],' jmin=', JMIN_IN[lp], JMAX_IS[lp]

      iimax=100L
      if lpPlot eq 120L then iimax=51L
      xtmp=fltarr(iimax)
      xtmp1=fltarr(iimax)
      xtmp2=fltarr(iimax)
      ytmp=fltarr(iimax)
;SH
      for ii=is,is-iimax+1,-1  do begin
         itmp=is-ii
         print,FORMAT='(I6,i6,e12.4,2f8.1,2f8.2)',ii,itmp,on_m3[ii,mp],tn_k[ii,mp],z_km[ii],glat_deg[ii,mp],glon_deg[ii,mp]
;  plot, on_m3[*,mp],z_km[*]
         xtmp[itmp]=on_m3[ii,mp] *1.0E-6 ;m-3-->cm-3
         xtmp1[itmp]=n2n_m3[ii,mp] *1.0E-6 ;m-3-->cm-3
         xtmp2[itmp]=o2n_m3[ii,mp] *1.0E-6 ;m-3-->cm-3
         ytmp[itmp]=z_km[ii]
      endfor                    ;ii    
      
      DEVICE, RETAIN=2, DECOMPOSED=0
      WINDOW,2,XSIZE=1000,YSIZE=800
      !p.multi=[0,1,1,0]
      

;plot vertical profile
      plot,xtmp,ytmp,/xlog $
           ,yrange=[90., YrangeMax] , xstyle=1$
           ,Xrange=[XrangeMin,1.E11], ystyle=1 $
                   ,charsize=2. $
                   ,title='20170413: '+whichModel+' O: glat='+STRTRIM( string(glat_deg[is,mp], FORMAT='(f7.0)'), 1)+' glon='+STRTRIM( string(glon_deg[is,mp], FORMAT='(f7.0)'), 1)+' LT='+STRTRIM( string(ltmp, FORMAT='(f6.1)'), 1)
           
           xyouts, 0.65, 0.9 $
                   ,'O'  $
                   , charsize=2.0, charthick=2.0, /norm, /noclip
           
           loadct,33
;SH [N2]
           oplot,xtmp1,ytmp $
                 ,linestyle=0 $
                 ,color=250     ;red
;SH [o2]
           oplot,xtmp2,ytmp $
                 ,linestyle=0 $
                 ,color=50      ;blue
           

;NH
           for ii=in,in+iimax-1,+1  do begin
              itmp=ii-in
              print,FORMAT='(I6,i6,e12.4,2f8.1,2f8.2)',ii,itmp,on_m3[ii,mp],tn_k[ii,mp],z_km[ii],glat_deg[ii,mp],glon_deg[ii,mp]
;  plot, on_m3[*,mp],z_km[*]
              xtmp[itmp]=on_m3[ii,mp]*1.0E-6 ;m-3-->cm-3
              xtmp1[itmp]=n2n_m3[ii,mp]*1.0E-6 ;m-3-->cm-3
              xtmp2[itmp]=o2n_m3[ii,mp]*1.0E-6 ;m-3-->cm-3
              ytmp[itmp]=z_km[ii]
           endfor               ;ii
;NH plot vertical profile [o]
           oplot,xtmp,ytmp $
                 ,linestyle=5
           


           loadct,39
;NH [N2]
           print,'min N2',min(xtmp1)
           oplot,xtmp1,ytmp $
                 ,linestyle=5 $
                 ,color=250     ;red

;NH [o2]
           print,'min o2',min(xtmp2)
           oplot,xtmp2,ytmp $
                 ,linestyle=5 $
                 ,color=50      ;blue
           

;ref1 alt=120km
           ytmp1=findgen(iimax)*0.+120.
           oplot,xtmp2,ytmp1 $
                 ,linestyle=1

;ref2 alt=470km
           ytmp1=findgen(iimax)*0.+470.
           oplot,xtmp2,ytmp1 $
                 ,linestyle=1

           xyouts, 0.45, 0.85 $
                   ,'N2'  $
                   , charsize=2.0, charthick=2.0, color=250, /norm, /noclip

           xyouts, 0.25, 0.7 $
                   ,'O2'  $
                   , charsize=2.0, charthick=2.0, color=50, /norm, /noclip


           output_png, whichModel+'profiles'+'_mlat'+STRTRIM( string(mlat_deg[in], FORMAT='(f6.0)'), 1)+'.20170413.png'
           STOP
endif 
endif;sw_debug_prfl eq 1 then begin


       if ( sw_frame eq 0 ) then begin ;magnetic
          plot_yy[mp,lps] = mlat_deg[i]
          plot_xx[mp,lps] = mlon_deg[mp]
          if (xmax360 eq 0) AND (mlon_deg[mp] ge +180.) then plot_xx[mp,lps]=mlon_deg[mp]-360.
       endif else if ( sw_frame eq 1 ) then begin ;geographic
          plot_yy[mp,lps] = glat_deg[i,mp]
          plot_xx[mp,lps] = glon_deg[i,mp]
          if (xmax360 eq 0) AND (glon_deg[i,mp] ge +180.) then plot_xx[mp,lps]=glon_deg[i,mp]-360.

;       endif else if ( sw_frame eq 2 ) then begin ;noon at center
;          plot_yy[mp,lps] = mlat_deg[i];
;
;;dbg20140923: 129 is jicamarca
;  ltime = ut_hr + glon_deg[i,mp]/15.
;  if ( ltime ge 24. ) then ltime = ( ltime MOD 24. )
;;extract near noon data
;;if ( lp eq 129 ) then  begin
;;  if ltime ge ltimemin and ltime le ltimemax then  $
;;    print, 'SH: mp',mp,'lp',lp,'LT',ltime, ut_hr, mlat_deg[i],glon_deg[i,mp]
;;  mp_output_noonprfl = mp
;;  glon_output = glon_deg[i,mp]
;;endif
;
;          plot_xx[mp,lps] = ltime

       endif ;( sw_frame eq 0 ) then begin ;magnetic


;20131209output to file:(2)SH
;20140808ipe 1st paper
; bug??? output for SH only??? why???
       if ( sw_output2file_ascii eq 1 ) AND ( ABS( mlat_deg[i]) le mlatmax ) then begin
          ltime=UT_hr + glon_deg[i,mp]/15.
          if ( ltime ge 24. ) then ltime = ( ltime MOD 24. )
          if ( ltime ge ltimemin) AND ( ltime le ltimemax ) then begin
             ncount=ncount+1
             printf, luntmp, UT_hr, plot_zz[mp,lps], TE_TI_k[3-1,i,mp], TE_TI_k[1-1,i,mp],mlat_deg[i], glon_deg[i,mp], ltime 
          endif                 ;( ltime ge 11.) AND ( ltime le 15. ) then begin
       endif                    ;( sw_output2file_ascii eq 1 ) then 


       BREAK; exit from the i loop
    endif; ( z_km[i] le ht_plot ) AND ( z_km[i+istep] gt ht_plot )  then begin 

;nm20130523: SH added nmf2/hmf2
   endif else begin ; ( VarType ge 7 ) then begin
     for jth=0,ISPEC-1 do  nel[i]=nel[i]+XIONN_m3[jth,i,mp]

     if ( i eq midpoint ) then  begin
       result = MAX( nel,Max_Subscript ) ;find NmF2
;       print,'SH: NmF2', result,Max_Subscript, nel[Max_Subscript],Z_km[Max_Subscript],mp,lp
       if ( VarType eq 7 ) then $
         plot_zz[mp,lps] = nel[Max_Subscript] * factor $ ;NmF2
       else if ( VarType eq 8 ) then $
         plot_zz[mp,lps] = Z_km[Max_Subscript]    ;hmF2

       if ( sw_frame eq 0 ) then begin ;magnetic
          plot_yy[mp,lps] = mlat_deg[Max_Subscript]
          plot_xx[mp,lps] = mlon_deg[mp]
          if (xmax360 eq 0) AND (mlon_deg[mp] ge +180.) then plot_xx[mp,lps]=mlon_deg[mp]-360.
       endif else if ( sw_frame eq 1 ) then begin ;geographic
          plot_yy[mp,lps] = glat_deg[Max_Subscript,mp]
          plot_xx[mp,lps] = glon_deg[Max_Subscript,mp]
          if (xmax360 eq 0) AND (glon_deg[Max_Subscript,mp] ge +180.) then plot_xx[mp,lps]=glon_deg[Max_Subscript,mp]-360.
;dbg20140923 SH
       endif else if ( sw_frame eq 2 ) then begin ;noon in the middle
;dbg20140923: 129 is jicamarca
;if ( lp eq 129 ) then  begin
  ltime = ut_hr + glon_deg[Max_Subscript,mp]/15.
  if ( ltime ge 24. ) then ltime = ( ltime MOD 24. )
;  if ltime ge ltimemin and ltime le ltimemax then  $
;    print, 'SH: mp',mp,'lps',lps,'LT',ltime, ut_hr
;endif

         plot_yy[mp,lps] = mlat_deg[Max_Subscript]
         plot_xx[mp,lps] = ltime
;nm20150224         plot_zz1[mp,lps] = Vn_ms1[2-1,Max_Subscript,mp]*fac_wind




       endif ;( sw_frame eq 0 ) then begin ;magnetic


;dbg20140825 validate wind direction nmf2
;if ( mp eq 51 ) then print, mp,lp, lps
if ( sw_output_nmf2 eq 2 AND mp eq 51 AND lp eq 36 ) then begin
if ( n_read eq 0 ) then begin
;   luntmpS=100
   flnmtmpS='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r319/trunk/run/'+rundir+'/nmf2SH.dat'
   openw,luntmpS,flnmtmpS, /GET_LUN
   print, 'SH nmf2 file created:',flnmtmpS
endif
;dbg20140825
     ltime=UT_hr + glon_deg[i,mp]/15.
     if ( ltime ge 24. ) then ltime = ( ltime MOD 24. )
print, mp,lp,UT_hr, plot_zz[mp,lps], Z_km[Max_Subscript], glat_deg[Max_Subscript,mp], glon_deg[Max_Subscript,mp], ltime
Printf, luntmpS, UT_hr, plot_zz[mp,lps], Z_km[Max_Subscript], glat_deg[Max_Subscript,mp], glon_deg[Max_Subscript,mp], ltime
endif ;( lp eq 47 AND mp eq 49 ) then begin

  
    endif;     if ( i eq midpoint ) then 
  endelse ;( VarType lt 7 ) then begin





 endfor                         ; i=is,midpoint, istep  do begin

endfor ;lp=0,NLP-1 do begin

;if ( n_read eq 0 ) then begin
;   flnmtmpS='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r319/trunk/run/'+rundir+'/nmf2SH.dat'
;   openw,luntmpS,flnmtmpS, /GET_LUN
;   print, 'nmf2@12LT file created:',flnmtmpS
;endif
;i=midpoint ;idealy IS of the flux tube at 20mlat
;ltime=UT_hr + glon_deg[i,mp]/15.
;          if ( ltime ge 24. ) then ltime = ( ltime MOD 24. )
;          if ( ltime ge ltimemin) AND ( ltime le ltimemax ) then begin
;             ncount=ncount+1
;             printf, luntmp, UT_hr, ltime, glon_deg[i,mp]
;             print, UT_hr,ltime,glon_deg[i,mp] 
;             for ii=0,ny_max  DO printf, ii, plot_yy[mp,ii],plot_zz[mp,ii]

;RETURN


endfor ;mp=0,NMP-1 do begin

;nm20140926 output noon profile
;if ( sw_frame eq 2 AND sw_output_nmf2noon eq 1 ) then begin
;   if ( n_read eq 0 ) then begin
;;   luntmpN=100
;      flnmtmpN='/scratch1/portfolios/NCEPDEV/swpc/noscrub/Naomi.Maruyama/r319/;trunk/run/'+rundir+'/noon_nmf2.dat'
;      openw,luntmpN,flnmtmpN, /GET_LUN
;      print, 'noon nmf2/wind file created:',flnmtmpN
;   endif                        ;if n_read eq 0
;   printf, luntmpN, mp_output_noonprfl 
;print, 'check mp_output=', mp_output_noonprfl 
;   printf, luntmpN, ut_hr, glon_output 
;print, 'ut=', ut_hr,' glon=', glon_output 
;   printf, luntmpN, plot_xx[mp_output_noonprfl,*] ;lt
;print, 'ltime=', plot_xx[mp_output_noonprfl,129] ;lt
;   printf, luntmpN, plot_yy[mp_output_noonprfl,*] ;mlat
;   printf, luntmpN, plot_zz[mp_output_noonprfl,*] ;nmf2
;   printf, luntmpN, plot_zz1[mp_output_noonprfl,*] ;wind positive NORTHward
;endif                                              ;if sw_frame eq 2

;(2) when time = time_max 
; plotting
if ( sw_range eq 1 ) then begin
      if ( VarType eq 0 ) then begin 
;f107-180
;        zmin=-8.0E+10
        zmin=0.;1.e+10
;        zmax=7.5e+11
;        zmax=2.5e+1 ;E-region
;        zmax=2.e+12 ;F-region
        zmax=22.55 ;F-region
;        zmax=9.0E+11 ;dbg200km
;        zmax=5.58e+12
;92km
;        zmax=2.00e+11
;        zmin=1.00e+04
;f107-80
;        zmin=1.e+10
;        zmax=1.77e+12
      endif else if ( VarType ge 1 ) AND (VarType le 2 ) then begin 
        zmin=0.
        zmax=3000.
      endif else if ( VarType ge 3 ) AND ( VarType lt 4 ) then begin 
;        zmin=0.
;        zmax=8.
;        zmin=1.
;        zmax=2.5e+1
        zmin=0.0
        zmax=3.0;2.e+12 ;F-region
      endif else if ( VarType eq 4 ) then begin 
        zmin=200. ;tn
        zmax=1000.
      endif else if ( VarType eq 5 ) then begin 
        zmin=-300.;zonal +eastward
        zmax=+300.
      endif else if ( VarType eq 6 ) then begin 
;        zmin=+14.0 ;350km
;        zmax=+14.57;

        zmin=+17.;100km
        zmax=+18.;
      endif else if ( VarType eq 7 ) then begin 
        zmin=0.  ; 0.0084; 1.0e+10
        zmax= $
;40. ;for TW
22.55
;7.521 ;1.983e+12
      endif else if ( VarType eq 8 ) then begin 
        zmin=100.
        zmax=500.
      endif


;nm20150205 output for paper
;if ( sw_frame eq 2 AND sw_output_nmf2noon eq 1 ) then begin
;t   printf, luntmpN,ut_hr
;t   ncnt=0L
;t   for k=0,ny_max-1  do begin
;t      if ( plot_yy[80-1,k] ne 0.0 AND plot_zz[80-1,k] ne 0.0 ) then begin
;t         ncnt = ncnt + 1
;t         printf, luntmpN, plot_yy[80-1,k], plot_zz[80-1,k]
;t;    print,  ncnt,         plot_yy[80-1,k], plot_zz[80-1,k]
;t      endif
;t   endfor
;t   print, 'ncnt ', ncnt 
;endif                           ;if ( sw_output2file_ascii eq 1 ) then begin


;if ( sw_debug eq 1 ) then $  
      maxZ=MAX(plot_zz)
      minZ=MIN(plot_zz)

;dbg20170427

j=100L
jj=170+j
for i=60,70 do begin
  print,'NH',i,j,  plot_zz[i,j],plot_xx[i,j],plot_yy[i,j]
  print,'SH',i,jj, plot_zz[i,jj],plot_xx[i,jj],plot_yy[i,jj]
endfor


	where_result=where ( plot_zz gt zmax, count ) 
;	print, '1. where_result', where_result,'count',count
	for i=0,count-1 do $
	  plot_zz[ where_result[i] ]=zmax
;  	where_result=where ( plot_zz gt zmax, count ) 
;	print, '2. where_result', where_result,'count',count

endif else if ( sw_range eq 0 ) then begin
  zmax = max(plot_zz)
  zmin = min(plot_zz)
  maxZ=zmax
  minZ=zmin
endif
print,'maxZ=',MAXz,' minZ=',MINz

MAX_xymin=zmax ;1.0E+3
n_levels=100L

if ( xmax360 eq 1 ) then begin
  X_max=+360.
  X_min=+  0.
endif else if ( xmax360 eq 0 ) then begin
  X_max= +180.
  X_min= -X_max
endif
if ( sw_frame eq 2 ) then begin
x_max=24.
x_min=0.
endif

Y_max= +80.0
Y_min= -Y_max
if ( sw_frame eq 0 ) or ( sw_frame eq 2 ) then $
  MAG_GEO='magnetic' $
else if ( sw_frame eq 1 ) then $
  MAG_GEO='geographic'

if ( n_read eq 0 ) then begin 
  X_TITLE=MAG_GEO+' longitude[deg]'
  Y_TITLE=MAG_GEO+' latitude[deg]'

  if ( sw_frame eq 2 ) then   X_TITLE='local time [hr]'

endif else begin
  X_TITLE=' '
  Y_TITLE=' '
endelse

if ( sw_plot_contour eq 1 ) then begin
	text_color=255.
char_size=1.0
char_thick=1.0
n_ldct=33 ;39;5;39
;if ( vartype eq 6 ) then N_LDCT=70

	iwindow=1L
if  ( n_read eq 0 ) then begin 
	DEVICE, RETAIN=2, DECOMPOSED=0
	WINDOW,iwindow,XSIZE=1100*fac_window,YSIZE=1000*fac_window
;                   columns,rows
nColumnsMulti = getenv('nColumnsMulti')  ;for quick plot
nRowsMulti    = getenv('nRowsMulti'   )  ;for quick plot
	!p.multi=[0,nColumnsMulti,nRowsMulti,0] ;
;	!p.multi=[0,4,7,0] ;6dy
;	!p.multi=[0,4,7,0] ;3dy
;	!p.multi=[0,6,5,0] ;1dy wam-ipe
;	!p.multi=[0,4,5,0]
;!	!p.multi=[0,2,2,0] ;1hr

	loadct,n_ldct
endif  ;( n_read eq 0 ) then begin 
endif ;( sw_plot_contour eq 1 ) then begin

X0=0.10
X1=0.90
Y0=0.10
Y1=0.79
if ( sw_plot_grid ge 1 ) then begin
loadct,0

if ( sw_plot_grid eq 1 ) then begin
plot,plot_xx,plot_yy $
,xrange=[X_min,X_max], xstyle=1 $
,yrange=[Y_min,Y_max], ystyle=1  $
,XTITLE=X_TITLE,YTITLE=Y_TITLE $
,TITLE=VarTitle[VarType]+unit[VarType]+'  ht='+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+'km  UT[hr]='+STRTRIM( string(ut_hr, FORMAT='(F6.2)'),1 )+'.'+TEST $
,POSITION=[X0,Y0,X1,Y1] $
,COLOR=text_color $
,charsize=char_size,charthick=char_thick $
,PSYM=3, SYMSIZE=5.5 ;$
title_plr='ret'

endif  else if ( sw_plot_grid eq 2 ) then begin
;      plot_yy[mp,lp] = mlat_deg[i]
LPmax=30
rmax=90.- plot_yy[0,LPmax-1]
print, 'rmax',rmax
r=fltarr(NMP,LPmax)
theta=fltarr(NMP,LPmax)
for mp=0,NMP-1 do begin
for lp=0,LPmax-1 do begin
r[mp,lp] = 90.- plot_yy[mp,lp] ;degrees from the North mag pole
theta[mp,lp] =  plot_xx[mp,lp] * !PI /180.
if ( r[mp,lp] gt rmax ) then  begin
  r[mp,lp]=0.
  theta[mp,lp]=0.
endif
endfor
endfor
plot,/POLAR, r, theta $
,XSTY=4, YSTY=4 $
;,xrange=[X_min,X_max], xstyle=1 $
;,yrange=[Y_min,Y_max], ystyle=1  $
;,XTITLE=X_TITLE,YTITLE=Y_TITLE $
,TITLE='POLAR PLOT in MAGNETIC COORDINATE' $ ;VarTitle[VarType]+unit+'  ht='+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+'km  UT[hr]='+STRTRIM( string(ut_hr, FORMAT='(F6.2)'),1 ) $
,POSITION=[X0,Y0,X1,Y1] $
,COLOR=text_color $
;,charsize=char_size,charthick=char_thick $
,PSYM=3, SYMSIZE=5.5 ;$
AXIS, 0,0,XAX=0 $
,Xrange=[-rmax,+rmax], Xstyle=1
AXIS, 0,0,YAX=0 $
,Yrange=[-rmax,+rmax], Ystyle=1
title_plr='plr'

endif  ;else if ( sw_plot_grid eq 2 ) then begin

if ( sw_output2file eq 1 ) then begin

   if ( sw_frame eq 0 ) then $  ;magnetic
      title_frame='mag' $
   else if ( sw_frame eq 1 ) then $ ;geographic
      title_frame='geo' $
   else if ( sw_frame eq 2 ) then $ ;noon in  the center
      title_frame='lt'

Filename_png=plot_DIR+'ipe_grid_ht'+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+title_plr+'.'+title_frame+'.png'
print,'output to file=',filename_png
output_png, Filename_png
endif

LOADCT, n_ldct
STOP
RETURN
endif ;( sw_plot_grid eq 1 ) then begin

;debug
;mp=0
;print,'mp==',mp
;for l=0,nlp-1 do print,plot_zz(mp,l),plot_xx(mp,l),plot_yy(mp,l)




if ( sw_plot_contour eq 1 ) then begin
ut_hr_disp= ut_hr MOD 24.


LOADCT, 0
contour,plot_zz,plot_xx,plot_yy $
,/irregular $
,/fill $
,levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
,xrange=[X_min,X_max], /xstyle $
,yrange=[Y_min,Y_max], /ystyle $
,XTITLE=X_TITLE,YTITLE=Y_TITLE $
;,TITLE=VarTitle[VarType]+unit[VarType]+'  ht='+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+'km  UT[hr]='+STRTRIM( string(ut_hr, FORMAT='(F6.2)'),1 )+'_'+TEST $
,TITLE='UT '+STRTRIM( string(ut_hr_disp, FORMAT='(F8.4)'),1 )+'  MAX='+STRTRIM( string(maxZ, FORMAT='(F8.2)'),1 )+' MIN='+STRTRIM( string(minZ, FORMAT='(F8.2)'),1 ) $
;,POSITION=[X0,Y0,X1,Y1] $
,COLOR=text_color $
,charsize=char_size,charthick=char_thick $
,MAX_VALUE= MAX_xymin $
,/NODATA

LOADCT, n_ldct
contour,plot_zz,plot_xx,plot_yy $
,/irregular $
,/fill $
,levels=findgen(n_levels)*(zmax-zmin)/float(n_levels-1) +zmin $
,xrange=[X_min,X_max], /xstyle $
,yrange=[Y_min,Y_max], /ystyle $
,XTITLE=X_TITLE,YTITLE=Y_TITLE $
;,TITLE=VarTitle[VarType]+unit[VarType]+'  ht='+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+'km  UT[hr]='+STRTRIM( string(ut_hr, FORMAT='(F6.2)'),1 )+'_'+TEST $
,TITLE='UT '+STRTRIM( string(ut_hr_disp, FORMAT='(F6.2)'),1 ) $
;,POSITION=[X0,Y0,X1,Y1] $
,COLOR=text_color $
,charsize=char_size,charthick=char_thick $
,MAX_VALUE= MAX_xymin $
,/OVERPLOT




if ( sw_debug eq 1 ) then  print,'MAX=',MAX(plot_zz),' MIN=',MIN(plot_zz)
; add MIN & MAX values
;xyouts, 0.5, 0.84 $
;, 'MIN='+STRTRIM(STRING( MIN(plot_zz), FORMAT='(E11.3)'),1)+' MAX='+STRTRIM(STRING( MAX(plot_zz), FORMAT='(E11.3)'),1)  $
;, charsize=1.0, charthick=1.0, /norm, /noclip

if ( n_read eq n_plt_min ) then begin
xyouts, 0.35, 0.04 $
,VarTitle[VarType]+unit[VarType]+'  ht='+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+'km'+' '+input_DIR0 $
, charsize=0.85, charthick=0.8, /norm, /noclip

charsize_colorbar=4.0
if VarType eq 6 then $
  format_colorbar='(f7.2)' $
else $
  format_colorbar='(E9.1)'

font=1 ;true-type 

dX=0.40
dY=0.010
X00=0.30
;y00=0.17;1day run
y00=0.08 ;3-6dy run
position=[X00, Y00, (X00+dX), (Y00+dY)] ;for horizontal bar
COLORBAR, BOTTOM=bottom, CHARSIZE=charsize_colorbar, COLOR=color, DIVISIONS=divisions $
        , FORMAT=format_colorbar, POSITION=position, MAXRANGE=zmax,MINRANGE=zmin $
        , NCOLORS=ncolors,TITLE=title,VERTICAL=vertical,TOP=top,RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
        , _EXTRA=extra, INVERTCOLORS=invertcolors, TICKNAMES=ticknames
endif ;( n_read eq n_plt_min


print,'n_read=',n_read,' n_plt_max=', n_plt_max
if ( n_read eq n_plt_max $
) AND ( sw_output2file eq 1 ) then begin
   if ( sw_frame eq 0 ) then $  ;magnetic
      title_frame='mag' $
   else if ( sw_frame eq 1 ) then $ ;geographic
      title_frame='geo'  $
   else if ( sw_frame eq 2 ) then $ ;noon in the center
      title_frame='lt'
Filename_png= $
plot_DIR+rundir+'_'+VarTitle[VarType]+'_ht'+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+title_frame+'.quick.png'
;plot_DIR+'quick/'+TEST+'_'+rundir+'_'+VarTitle[VarType]+'_ht'+STRTRIM( string(ht_plot, FORMAT='(F4.0)'),1 )+title_frame+'.quick.png'
output_png, Filename_png
endif ;( sw_output2file eq 1 ) then begin

endif ;( sw_plot_contour eq 1 ) then begin


endfor ;VarType=1,1 do begin
end ;pro contour_lon_lat
