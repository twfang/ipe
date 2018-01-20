;density profile along the ERG orbit 
;!!!CAUTION! IPE density always at equatorial plane
;v4: L & MLT factors included in interpolation
pro restore_erg_v4

sw_debug=0L
titlePlot2='22.8-24UT ERG MLT'
swPlotProfile=1
filenameSav='/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/fig/obana/erg2017-09-08.sav'
restore,filename=filenameSav
print,'MIN UT=',MIN(utSav),'MAX UT=',MAX(utSav)
;---read ERG orbit
iReadMax=190L
uTime=fltarr(iReadMax)
Rre=fltarr(iReadMax)
glat=fltarr(iReadMax)
mlt=fltarr(iReadMax)
Lval=fltarr(iReadMax)
Bt=fltarr(iReadMax)
rPlot=fltarr(iReadMax)
zPlot=fltarr(iReadMax)

openr,lun3000,'/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/mpi20160330v2/run17/erg/ERG_orbit04pos.txt',/get_lun
iDate0=0L
tmpString='# Date    Time[H]   R[RE]    GLAT    MLT[H] L-value   B[nT]'
iRead=-1L
   readf,lun3000,tmpString
   print,'line#1=',tmpString
   readf,lun3000,tmpString
   print,'line#2=',tmpString
   while ( eof(lun3000) eq 0 ) do begin
      readf,lun3000,iDate0,uTime0,rre0,glat0,mlt0,lval0,bt0,format='(i8,f10.6,f9.6,f8.3,f7.3,f10.6,f8.1)'
      iRead=iRead+1
      print,iRead,iDate0,uTime0,rre0,glat0,mlt0,lval0,bt0,format='(i4,i9,f10.6,f9.6,f8.3,f7.3,f10.6,f8.1)'
      uTime[iRead]=uTime0
      Rre[iRead]=rre0
      glat[iRead]=glat0
      mlt[iRead]=mlt0
      Lval[iRead]=lval0
      Bt[iRead]=bt0

   endwhile ;eof(lun3000
   print,'iReadMax=',iRead


   Re=6.3712E+03                ;km                                                         
   r_ref=Re+90.                 ;km for reference ht
;---search for matching UT & MLT

;   for i=0,iReadMax-1 do begin
   for i=0,63 do begin ;until 9/8 23.99UT
      print,'i=',i,' uTime=',uTime[i],' mlt=',mlt[i]

;      for j=0,nMax-1 do begin
      for j=0,95 do begin
         print,'j=',j,' utSav[j]=',utSav[j]
         if ( utSav[j] le uTime[i] and utSav[j+1] gt uTime[i] ) then begin
            print,'UT MATCH! ',j,utSav[j],uTime[i],utSav[j+1]


            for m=0,nmp-1 do begin

                 mp0=m
                 mp1=m+1

                 if mp1 ge nmp then  mp1=mp1-nmp

if sw_debug eq 1 then              print,'m=',mp0,' mltSav[j,m]=',mltSav[j,mp0]
              if ( mltSav[j,mp0] le mlt[i] and mltSav[j,mp1] gt mlt[i] ) then begin
                 print,mp0,'MLT MATCH! mltSav=',mltSav[j,mp0],mlt[i],mltSav[j,mp1]

                 mltFac = ( mlt[i]-mltSav[j,mp0] ) / ( mltSav[j,mp1]-mltSav[j,mp0] )


;search for matching L: 1.5<=L<=6
;lp=21 L=6.19546
                 for lp=21, nlp-2 do begin

;calculate L value at lp
;lp0
                    iin=jmin_in[lp]-1L
                    midpoint = JMIN_IN[lp]-1 + ( JMAX_IS[lp] - JMIN_IN[lp] )/2
                    sinthet = SIN( (90.- mlat_deg[iin])  *!PI / 180. )
                    lval0    = r_ref / ( Re * sinthet * sinthet )
if sw_debug eq 1 then   print,Lval[i],' Lval=',lval0,' lp0=',lp,' mlat=',mlat_deg[iin],' z_km=',z_km[iin]
;lp1
                    iin=jmin_in[lp+1]-1L
                    midpoint = JMIN_IN[lp+1]-1 + ( JMAX_IS[lp+1] - JMIN_IN[lp+1] )/2
                    sinthet = SIN((90.- mlat_deg[iin])  *!PI / 180. )
                    lval1    = r_ref / ( Re * sinthet * sinthet )

;stop

                    if ( lval1 le Lval[i] and lval0 gt Lval[i] ) then begin
                       lp0=lp
                       lp1=lp+1
                       lFac = ( Lval[i]-lval1 ) / ( lval0-lval1 )
;error trap
                       if lFac lt 0.0 or lFac gt 1.0 then begin 
                          print,'!STOP! INVALID lFac=',lFac
                          STOP
                       endif
                       break    ;from lp loop
                    endif       ;lval0
                 endfor         ;lp
;--after match is found:
; at the next level of interpolation, MLT factor should be included!
                 zPlot[i] = zSav[j,m,lp1] + ( zSav[j,m,lp0] - zSav[j,m,lp1] ) * lFac
                 rPlot[i] = rSav[j,m,lp1] + ( rSav[j,m,lp0] - rSav[j,m,lp1] ) * lFac

                 

              endif;mltsav[ 
          
           endfor               ;m=0,nmp-1 do begin
         endif                  ;utSav[j] le
      endfor                    ;j=0,nMax do begin



JUMP1: print,'end i loop:i=',i,uTime[i]
endfor ;i=0,iReadMax do begin

;---plotting
;if swPlotProfile eq 1 then begin
      m2cm=1.0e+02
      xmax=6.
      xmin=1.5
      col_max = 255.9999999999999999
   axis_color =255.
   char_size=1.

;nReadMin=0L
;for nRead=nReadMin,nMax-1 do begin
   ;for mp=0,nmp-1 do begin

      ;if ( i eq 0 ) then begin

         iwindow=2L
         DEVICE, RETAIN=2, DECOMPOSED=0
         WINDOW,0,XSIZE=700,YSIZE=500

         n_ldct=0               ;black+white                                                         
         loadct, n_ldct

         plot,rPlot,zPlot $
              ,xrange=[xmin    ,xmax    ], xstyle=1 $
              ,yrange=[0.5e+00,1.e+04], ystyle=1 $
              ,/YLOG $
              ,title=titlePlot+' '+titlePlot2+' Ne number density: ERG v4'  $
              ,linestyle = 0 $
              ,color=axis_color $
              ,charsize=char_size $
              ,/NODATA
      ;endif; (nRead eq nReadMin AND mp eq 0 ) then begin


if sw_debug eq 1 then for i=0,63 do print, i,uTime[i],rplot[i],zplot[i]

     ; n_ldct=39                  ;black+white                                                        
     ; loadct, n_ldct
      oplot,rPlot,zPlot $
            ,linestyle=0 $
            ,color = 255


;MLT symbol
;      n_ldct=39                  ;black+white                                                        
;      loadct, n_ldct
;for i=0,63 do $
;      oplot,rPlot[i],zPlot[i] $
;;            ,linestyle=0 $
;      ,psym=3,symsize=3 $
;            ,color = mlt[i]*col_max/24.


  ; endfor                       ; mp=0,nmp-1 do begin
;endfor   ; nRead=nReadMin,nMax-1 do begin

;endif                           ;swPlotProfile eq 1 then begin




print,'pro restore_erg_v4 finished!'
end                             ;pro restore_erg_v4
