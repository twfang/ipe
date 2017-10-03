pro rd_pe

  mplt=1L
  lplt=21L ;1L
  machine= $
;'mac'
'theia'
sw_output2file=1L
swLun17=1L
sw_debug=0L
fac_window=0.8

if machine eq 'theia' then begin
   ppath='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/runs/tmp20150730/trunk/run/'
   rundir= $
;'1481414104_ipe_theia_intel_parallel2_80' ;v30 new run without perp trans
'1481593483_ipe_theia_intel_parallel2_80';v31 new run without conjugate PE
   pltDir='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/fig/20161128kitamura/'
endif else if machine eq 'mac' then begin
   ppath='.'
   rundir=''
endif

input_flnm='fort.166'
input_flnm1='fort.167'
input_flnm2='plasma17' ;sza
lun166=166L
lun167=167L
lun17=17L
flnm=ppath+rundir+'/'+input_flnm
print,'166opening file',flnm
openr, lun166, flnm, /get_lun
flnm=ppath+rundir+'/'+input_flnm1
print,'167opening file',flnm
openr, lun167, flnm, /get_lun
flnm=ppath+rundir+'/'+input_flnm2
print,'17opening file',flnm
if swLun17 eq 1 then  openr, lun17, flnm, /get_lun  , /F77_UNFORMATTED
npts2d=31287L
nmp=80L
iemax=98L
sza=fltarr(npts2d,nmp)
energy=fltarr(iemax)
phiup=fltarr(iemax)
phidwn=fltarr(iemax)
;Conjugate
phiupC=fltarr(iemax)
phidwnC=fltarr(iemax)

format166='(i6,2i3)'
format1661='(F6.1,12E9.2)'
string_tmp166='E   SIGEXO  SIGIONO  SIGEXN2  SIGION2  SIGEL    PRED     FYSUM    TSIGNE   PHIUP    PHIDWN   PRODUP   PRODWN'
ut0=864000L
utime=0L
mp0=0L
lp0=0L
mp1=0L
lp1=0L

nmax=13L                        ;output every 2hr
jmax=22L                     ;=LpMAX

totalP=fltarr(nmax,jmax)
eav=fltarr(nmax,jmax)
sza_deg=fltarr(nmax,jmax)
sza_degs=fltarr(nmax,jmax)
for nt=0,nmax-1 do begin
; read loop
   for j=0,jmax-1 do begin

;loop energy
      readf, LUN166, utime,mp0,lp0, FORMAT=format166
      print,j,'166utime=',utime,' mp0=',mp0,' lp0=',lp0
      readf, LUN167, utime,mp1,lp1, FORMAT=format166
      print,j,'167utime=',utime,' mp1=',mp1,' lp1=',lp1
   



;energy loop energy
      for i=0,iemax-1 do begin
         ie=iemax-1L-i
         if sw_debug eq 1 then print,' i=',i, ' ie=', ie
         
;E(IE),SIGEX(1),SIGION(1),SIGEX(3),SIGION(3)
;     >    ,SIGEL(1),PRED(IWR),FYSUM(IWR),TSIGNE(IWR),PHIUP(IWR)
;     >    ,PHIDWN(IWR),PRODUP(IE,IWR),PRODWN(IE,IWR)
;fort166;NH
         readf,lun166, tmp0,tmp1,tmp2,tmp3,tmp4 $
               , tmp5,tmp6,tmp7,tmp8,tmp9 $
               , tmp10,tmp11,tmp12, FORMAT=format1661
         
         energy[ie]=tmp0   
         phiup[ie]=tmp9   
         phidwn[ie]=tmp10   
      
         if $
sw_debug eq 1 and $
            i eq 0 then print, tmp0,tmp1,tmp2,tmp3,tmp4 $
                               , tmp5,tmp6,tmp7,tmp8,tmp9 $
                               , tmp10,tmp11,tmp12 $
                               ,FORMAT=format1661
         
;fort167;SH conjugate
         readf,lun167, tmp0,tmp1,tmp2,tmp3,tmp4 $
               , tmp5,tmp6,tmp7,tmp8,tmp9 $
               , tmp10,tmp11,tmp12, FORMAT=format1661
      
         phiupC[ie]=tmp9   
         phidwnC[ie]=tmp10   
      
         if $
           sw_debug eq 1 and $
            i eq 0 then print, tmp0,tmp1,tmp2,tmp3,tmp4 $
                               , tmp5,tmp6,tmp7,tmp8,tmp9 $
                               , tmp10,tmp11,tmp12 $
                               ,FORMAT=format1661

      endfor                    ;i  energy

;evaluate surface area
;for i=0,iemax-1 do print, i, energy[i]
totalP[nt,lp0-1]=0.0E0
totalE=0.0E0
      for i=0,iemax-1 do begin
         if i eq 0 then $
            difE=energy[i]-0. $
         else $
            difE=energy[i]-energy[i-1]

         if sw_debug eq 1 then print,i,' dif=',difE
         difEPhi=difE*phiup[i]
         EPhi   =energy[i]*phiup[i]
         totalP[nt,lp0-1]=totalP[nt,lp0-1]+difEPhi 
         totalE=totalE+EPhi 
      endfor
      print,'photoelectron density flux [cm-2 s-1]=', totalP[nt,lp0-1]
      eav[nt,lp0-1]=totalE/totalP[nt,lp0-1]
      print,'average energy [eV]=', eav[nt,lp0-1]
         
      if swLun17 eq 1 AND j eq 0 then  readu,lun17, sza 

;read grid for IN/IS
      if nt eq 0 and j eq 0 then begin
         filename_grid_sav='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/grid/plt/plasma_grid.2xdyn.sav'
         print,'restoring the grid file=',filename_grid_sav
         restore, filename=filename_grid_sav
         if sw_debug eq 1 then print, jmin_in[lp0-1],jmax_is[lp0-1]
      endif

      iz=86L; i at 800km
;NH
      isza=iz + jmin_IN[lp0-1]-1
      sza_deg[nt,lp0-1]=(sza(isza,mp0)*180./!PI)
;SH
      iszaS=jmax_IS[lp0-1]-iz-1
      sza_degS[nt,lp0-1]=(sza(iszaS,mp0)*180./!PI)
      if sw_debug eq 1 then print,z_km[isza],z_km[iszaS],mp0,lp0
      print,'sza[deg]',sza_deg[nt,lp0-1],sza_degS[nt,lp0-1]


;start plotting
      if sw_debug eq 1 then $
         print,'check mplt', mplt,mp0,lplt,lp0
      if mplt eq mp0 AND lplt eq lp0 then begin
         x_min=0.
         x_max=80.
         y_min=0.
         y_max=8.0E+07
         text_color=255.
         char_size=1.5
         char_thick=1.5
         n_ldct=39

         iwindow=0
         DEVICE, RETAIN=2, DECOMPOSED=0
         WINDOW,iwindow,XSIZE=800*fac_window,YSIZE=800*fac_window
         !p.multi=[0,1,1,0]
         loadct,0
      

;plot photoelectron enerfy spectrum
         plot, energy, phiup, linestyle=0 $
            ,xrange=[X_min,X_max], /xstyle $
            ,yrange=[y_min,y_max], /ystyle $
            ,color=lineColor $
            ,title='photoelectron energy spectra at 800 km' $
            ,xtitle='photoelectron energy [eV]' $
            ,ytitle='photoelectron flux [eV-1 cm-2 sec-1]' $
            ,charsize=char_size,charthick=char_thick $
            ,/NODATA


      loadct,n_ldct
      line_style=0
      line_thick=4.
      line_color=250            ;red
      oplot, energy, phiup $
             , color=line_color $
             , linestyle=line_style $
             , thick=line_thick
      line_style=2
      oplot, energy, phidwn $
             , color=line_color $
             , linestyle=line_style $
             , thick=line_thick
      
;conjugate
      line_style=0
      line_thick=2.
      line_color=50             ;blue
      oplot, energy, phiupC $
             , color=line_color $
             , linestyle=line_style $
             , thick=line_thick
      line_style=2
      oplot, energy, phidwnC $
             , color=line_color $
             , linestyle=line_style $
             , thick=line_thick



      xyouts, 0.73, 0.009 $
              ,'ut='+STRTRIM( string((utime-ut0)/3600., FORMAT='(F6.0)'),1 )+' sza='+STRTRIM( string(sza_deg[nt,lp0-1], FORMAT='(F6.0)'),1 )+'; '+STRTRIM( string(sza_degS[nt,lp0-1], FORMAT='(F6.0)'),1 )+' Eav='+STRTRIM( string(eav[nt,lp0-1],  FORMAT='(f6.1)') ,1) $
              ,charsize=1., charthick=1., /norm, /noclip
      
      xyouts, 0.02, 0.009 $
              ,'mp='+STRTRIM( string(mp0, FORMAT='(i2)'),1 )+' lp='+STRTRIM( string(lp0, FORMAT='(i2)'),1 ) $
              ,charsize=1., charthick=1., /norm, /noclip
      
;Fig2 of Sun et al 1998
      filename_image=pltDir+rundir+'/Fig2photoelectronEnergySpectra'+'_ut'+STRTRIM( string((utime-ut0)/3600., FORMAT='(F6.0)'),1 )+'_mp'+STRTRIM( string(mp0, FORMAT='(i2)'),1 )+'lp'+STRTRIM( string(lp0, FORMAT='(i2)'),1 )+'.png'
      if sw_output2file eq 1 then begin
         print, filename_image
         output_png, filename_image
      endif
   endif                        ;mplt eq mp0 AND lplt eq lp0 then begin
      
   endfor                       ;jmax
endfor                          ;nt=0,nmax-1 do begin

free_lun,lun166
free_lun,lun167
if swLun17 eq 1 then   free_lun,lun17

; plot fig4
         iwindow=1
         DEVICE, RETAIN=2, DECOMPOSED=0
         WINDOW,iwindow,XSIZE=800*fac_window,YSIZE=800*fac_window
         !p.multi=[0,1,2,0]
for nt=0,nmax-1 do begin
   if nt eq 0 then   plot, sza_deg[nt,*], totalP[nt,*],/YLOG  $
            ,xrange=[0.,180.], /xstyle $
            ,yrange=[0.9E+5,4.E+9], /ystyle $
   else             oplot, sza_deg[nt,*], totalP[nt,*] 
endfor
for nt=0,nmax-1 do begin
   if nt eq 0 then   plot, sza_deg[nt,*], eav[nt,*]  $
            ,xrange=[0.,180.], /xstyle $
            ,yrange=[0.,20.], /ystyle $
   else             oplot, sza_deg[nt,*], eav[nt,*]  
endfor
;Fig4 of Su et al 1998
filename_image=pltDir+'Fig4_mp'+STRTRIM( string(mplt, FORMAT='(i2)'),1 )+'lp'+STRTRIM( string(lplt, FORMAT='(i2)'),1 )+'.png'
         output_png, filename_image



print, 'pro rd_pe finished'
end                             ;pro rd_pe
