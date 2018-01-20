;L profiles at a given MLT
pro restore_erg_v0

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

   endwhile
   print,'iReadMax=',iRead
;---search for matching UT & MLT
   rPlot=fltarr(nlp)
   zPlot=fltarr(nlp)
   for i=0,iReadMax-1 do begin
      print,'i=',i,' uTime=',uTime[i],' mlt=',mlt[i]

      for j=0,nMax-1 do begin
         print,'j=',j,' utSav[j]=',utSav[j]
         if ( utSav[j] le uTime[i] and utSav[j+1] gt uTime[i] ) then begin
            print,'UT MATCH! ',j,utSav[j],uTime[i],utSav[j+1]


            for m=0,nmp-1 do begin
              print,'m=',m,' mltSav[j,m]=',mltSav[j,m]
              if ( mltSav[j,m] le mlt[i] and mltSav[j,m+1] gt mlt[i] ) then begin
                 print,'MLT MATCH! ',m,mltSav[j,m],mlt[i],mltSav[j,m+1]
                 zPlot[*]=zSav[j,m,*]
                 rPlot[*]=rSav[j,m,*]
                 GOTO, JUMP1

              endif             ; 
          
           endfor               ;m=0,nmp-1 do begin
         endif                  ;utSav[j] le
      endfor                    ;j=0,nMax do begin



JUMP1: print,'start plotting:i=',i,uTime[i]
;---plotting
;if swPlotProfile eq 1 then begin
      m2cm=1.0e+02
      xmax=5.
      xmin=1.5
      col_max = 255.9999999999999999
   axis_color =255.
   char_size=1.

;nReadMin=0L
;for nRead=nReadMin,nMax-1 do begin
   ;for mp=0,nmp-1 do begin

      if ( i eq 0 ) then begin

         iwindow=2L
         DEVICE, RETAIN=2, DECOMPOSED=0
         WINDOW,0,XSIZE=700,YSIZE=500

         n_ldct=0               ;black+white                                                         
         loadct, n_ldct

         plot,rPlot[*],zPlot[*] $
              ,xrange=[xmin    ,xmax    ], xstyle=1 $
              ,yrange=[1.e+00,1.e+06], ystyle=1 $
              ,/YLOG $
              ,title=titlePlot+' '+titlePlot2+' Ne number density: ERG'  $
              ,linestyle = 0 $
              ,color=axis_color $
              ,charsize=char_size $
              ,/NODATA
      endif; (nRead eq nReadMin AND mp eq 0 ) then begin


      n_ldct=39                  ;black+white                                                        
      loadct, n_ldct
      oplot,rPlot[0:nlp-1],zPlot[0:nlp-1] $
            ,linestyle=0 $
            ,color = mlt[i] * col_max / 24.


  ; endfor                       ; mp=0,nmp-1 do begin
;endfor   ; nRead=nReadMin,nMax-1 do begin

;endif                           ;swPlotProfile eq 1 then begin
endfor ;i=0,iReadMax do begin

;---add colorbar
title='';colorbar title'                                                                             
font=1  ;True-Type: 1.                                                                               
; X_SIZE=29.7                                                                                         
; Y_SIZE=20.5                                                                                         
; X0 = 2.2                                                                                            
; Y0 = 2.6                                                                                            
; Y_cmargin=1.0                                                                                       
;dYc = 0.28                                                                                           
charsize_colorbar=2.
format_colorbar='(F4.0)'
;position=[X0/X_SIZE, (Y1+Y_cmargin)/Y_SIZE, X1/X_SIZE,
;(Y1+Y_cmargin+dYc)/Y_SIZE]                    
position=[0.15, 0.94, 0.940, 0.95]
zmin=0.0
zmax=24.0
COLORBAR, BOTTOM=bottom, CHARSIZE=charsize_colorbar, COLOR=color, DIVISIONS=divisions $
          , FORMAT=format_colorbar, POSITION=position, MAXRANGE=zmax, MINRANGE=zmin $
        , NCOLORS=ncolors, TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right $
        , MINOR=minor, RANGE=range, FONT=font, TICKLEN=ticklen $
          , _EXTRA=extra, INVERTCOLORS=invertcolors,  TICKNAMES=ticknames

;---

print,'pro restore_erg_v0 finished!'
end                             ;pro restore_erg_v0
