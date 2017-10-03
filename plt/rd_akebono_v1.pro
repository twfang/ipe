;date: 20170215
;purpose: restore template, read ascii file, and plot akebono data
;v0: 2013-03-16_ 10<ut<10.45
;v1: 2015-03-16
pro rd_akebono_v1
version='v2'
;2013-03-16 10ut
if version eq 'v0' then begin
cmpDate='20130316'
   iakMin=0L
   iakMax=13L
   x_min=10.+18./60.
   x_max=10.+25./60.
   y_min=3500.
   y_max=5000.
ModIak=0L

;2013-03-16 13ut   
endif else if version eq 'v1' then begin
cmpDate='20130316'
   iakMin=13L
   iakMax=30L
   x_min=14.+49./60.
   x_max=14.+58./60.
   y_min=3400.
   y_max=6025.
ModIak=1L

;2015-03-16 00ut   
endif else if version eq 'v2' then begin
cmpDate='20150316'
   iakMin=0L
   iakMax=17L
   x_min=0.+0./60.
   x_max=0.+9./60.
   y_min=8887.
   y_max=11150.
ModIak=0L
  
endif ;version eq 'v0' then begin

flnm_png='cmpAkebono'+cmpDate+'.'+version+'.png'
  sw_debug=0
  machine= $
    ; 'mac'
     'theia'
   ;'/home/naomi.maruyama/' ;naomilx
   ;

if machine eq 'theia' then $
   workdir='/home/Naomi.Maruyama/sandbox/ipe/trunk/plt/' $
else if machine eq 'mac' then $
   workdir='/Users/naomimaruyama/scrub/tmp/tmp20170215akebono/'

myfile=workdir+'ne-'+cmpDate+'.txt'
restore,workdir+'mytemplateAkebono.sav'
PLOT_ASCII = READ_ASCII(myfile, TEMPLATE = mytemplate)

;print, PLOT_ASCII 
if sw_debug eq 1 then  help, PLOT_ASCII 
if sw_debug eq 1 then  print, 'size=',size(plot_ascii)
  UTh=FIX(PLOT_ASCII.UTime*1.0e-4)
if sw_debug eq 1 then   print,'uth=', uth
  UTm= FIX(( PLOT_ASCII.UTime - UTh*1.0e+4 )*1.0e-2)
if sw_debug eq 1 then   print,'utmin=', utm
  UTs= FIX(PLOT_ASCII.UTime - (UTh*1.0e+4 +UTm*1.0e+2) )
if sw_debug eq 1 then   print,'utsec=', uts
  UTSec = UTs + UTm*60. + UTh*3600.
if sw_debug eq 1 then   print,'UTSec[sec]=', UTSec
if sw_debug eq 1 then print,'size_result utsec=',size(utsec)
size_result=size(utsec)
nmax=size_result[1]
print,'nmax=',nmax
  MLTak= PLOT_ASCII.GMLT
if sw_debug eq 1 then   print,'MLT[hr]=', MLTak
if sw_debug eq 1 then print,'size_result MLT=',size(MLTak)


  Mlat= PLOT_ASCII.GMLAT
if sw_debug eq 1 then   print,'MLAT[deg]=', MLat

  NeAk= PLOT_ASCII.NEPERCC
if sw_debug eq 1 then   print,'Ne[/cc]=', NeAk
;endfor

;read ipe grid
;read_grid, JMIN_IN,

if machine eq 'theia' then $
  restore, filename='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/grid/plt/plasma_grid.2xdyn.sav' $
else if machine eq 'mac' then $
  restore,filename='/Users/naomimaruyama/save/sandbox/ipe/plt.bk20150326/plasma_grid.2xdyn.sav'

if sw_debug eq 1 then begin
           print,'JMIN', JMIN_IN[0:1]
print,'JMAX', JMAX_IS[0:1]
print,'z_km', z_km[0],size(z_km)
print,'mlat_deg', mlat_deg[0],size(mlat_deg)
endif ;sw_debug


;IPE read_plasma
if machine eq 'theia' then $
  path='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/mpi20160330v2/run2/' $
else if machine eq 'mac' then $
  path=workdir

rpath=path+'1462618349_ipe_theia_intel_parallel2_93/'
openr,lun00,rpath+'plasma00',/get_lun, /F77_UNFORMATTED ;o+
openr,lun01,rpath+'plasma01',/get_lun, /F77_UNFORMATTED ;h+
openr,lun02,rpath+'plasma02',/get_lun, /F77_UNFORMATTED ;he+
openr,lun2013,rpath+'fort.2013',/get_lun ;, /F77_UNFORMATTED
openr,lun0,rpath+'ut_rec',/get_lun ;, /F77_UNFORMATTED

imax=iakMax-iakMin+1
xplt0=fltarr(imax);ut
yplt0=fltarr(imax);akebono ne
xplt1=fltarr(imax)
yplt1=fltarr(imax);ipe ne
yplt2=fltarr(imax) ;he+
yplt3=fltarr(imax) ;o+
zalt=fltarr(imax)
zmlat=fltarr(imax)
zmlt=fltarr(imax)
zmp=LONarr(imax)
zlp=LONarr(imax)


ncnt=0L
NPTS2D=31287L
NMP=80L
dum0=fltarr(NPTS2D,NMP) ;o+
dum1=fltarr(NPTS2D,NMP) ;h+
dum2=fltarr(NPTS2D,NMP) ;he+

iak=iakMin
while iak le iakMax-1  do begin
   print,'iak=',iak,' UTSec[iak]=',utSec[iak]



ut00=518400
ut=0L
recnumber=0L

;while ncnt le nMax  do begin
if (iak MOD 2 ) eq ModIak then begin
  print,'start reading plasma: iak=', iak
   while (eof(lun0) eq 0) do begin
      readf, lun0,recnumber, ut
      ncnt=ncnt+1
      utipe=(ut-ut00)*1.0
      if sw_debug eq 1 then print,ncnt,' rec#',recnumber,' ut=',ut,' utipe[s]=', utipe

;read sunlon
      readf, LUN2013,sunlons1
                                ;print,'sunlons1=',sunlons1

      readu, lun00,dum0 ;o+
      readu, lun01,dum1 ;h+
      readu, lun02,dum2 ;he+
      
      if utipe ge utSec[iak] then begin
         print, ncnt,'utipe=', utipe
         print, 'utSec=', utSec[iak],iak
         ;print,'sunlons1=',sunlons1
         break                  ;exit from the while loop
      endif                     ;utipe ge uth then begin

   endwhile ;eof
endif; (iak MOD 2 ) eq 0 then begin
   ;endwhile                     ;ncnt le nMax then begin
;print,'after endwhile: ncnt=',ncnt

;---calcualte mlt
mlt = fltarr(nmp)
ddeg=360./NMP
mlon_deg=findgen(nmp)*ddeg
if sw_debug eq 1 then print,'mlon_deg=', mlon_deg
for mp=0,nmp-1   do begin
   mlt[mp]    = mlon_deg[mp]/15.0D0 - sunlons1 * 12.0D0 / !PI   +12.0 ;[hr] !CORRECT!
  if ( mlt[mp] lt  0. ) then  mlt[mp] = mlt[mp] MOD 24.
  if ( mlt[mp] ge 24. ) then  mlt[mp] = mlt[mp] - 24.

;if ( i eq mp ) then 


  if mlt[mp] ge MLTak[iak] then begin
    mplt=mp
    zmp[iak-iakMin]=mp
    zmlt[iak-iakMin]=mlt[mp]
    break
  endif ;mlt[mp] ge MLTak[iak] then 

endfor
print,'mp=',(mplt+1),' mlt[hr]=', mlt[mplt], MLTak[iak]
;print, 'sunlons1=', sunlons1
;---
ipts1D=195  ;3002; 4002km
di=5
;---choose lp

;ipts2D=jmin_in[lp]+ipts1D
lp0=38L
dlp=3
dmlat=3.16 ;2.5
for lp=(lp0-dlp),(lp0+dlp)  do begin
  ;for i=(jmin_in[lp]+ipts1D-di),(jmin_in[lp]+ipts1D+di)  do begin
  for i=(jmax_is[lp]-ipts1D+di),(jmax_is[lp]-ipts1D-di),-1  do begin
    print,lp,(i-jmin_in[lp]),'z_km=', z_km[i],mlat_deg[i],dum1[i,mplt]
    if z_km[i] ge 4000. then begin
      ipts2D=i
      zalt[iak-iakMin]=z_km[i]
      break
    endif
  endfor ;i
  if ABS(mlat_deg[ipts2D]-mlat[iak]) le dmlat then  begin
    zmlat[iak-iakMin]=mlat_deg[ipts2D]
    lplt=lp
    zlp[iak-iakMin]=lp
    break
  endif
endfor ;lp
print, 'lplt=',lplt,'ipts2D', ipts2D,mlat[iak]

;---assign data to plotting arrays

  xplt0[iak-iakMin]=utSec[iak]/3600.
  yplt0[iak-iakMin]=NeAk[iak]

  xplt1[iak-iakMin]=xplt0[iak-iakMin];utipe/3600.
  yplt1[iak-iakMin]=(dum0[ipts2D,mplt]+dum1[ipts2D,mplt]+dum2[ipts2D,mplt])*1.0e-6 ;m-3-->cm-3
  yplt2[iak-iakMin]=(dum2[ipts2D,mplt])*1.0e-6 ;m-3-->cm-3  ;he+
  yplt3[iak-iakMin]=(dum0[ipts2D,mplt])*1.0e-6 ;m-3-->cm-3  ;o+

  iak=iak+1
endwhile                        ;iak le iakMax  do begin
print, 'iak=', iak

;---plot
  device, decomposed = 0 ,retain=2
  window, 0 ,XSIZE=800,YSIZE=800
  !P.BACKGROUND=255
;   columns, rows
!P.MULTI=[0,1,3]


;--
loadct, 0
Sym_Size=4
char_size=3.
char_thick=2.0
plot, xplt0, yplt0,color=0, psym=4,symsize=Sym_Size $
  ,xrange=[x_min,x_max], xstyle=1  $
  ,yrange=[y_min,y_max], ystyle=1  $
  ;,xtitle='UT [hr]' $
  ,ytitle='Number Density [cm-3]' $
  ,title='Comparison with Akebono:'+cmpDate $
, charsize = char_size, charthick = char_thick  $
  ,/NODATA

loadct, 39
oplot, xplt0, yplt0, color=50, psym=4,symsize=Sym_Size 
for i=(iakMin-iakMin),(iakMax-iakMin) do print,i, xplt1[i],yplt1[i],yplt2[i],yplt3[i],zmlt[i],zalt[i],zmlat[i] $
,zmp[i],zlp[i] $
,FORMAT='(i2,f7.2,3f8.1,f7.1,f7.0,f8.2,2i4)'
oplot, xplt1, yplt1, color=250, psym=1,symsize=Sym_Size

;2;mlat
plot, xplt0, zmlat,color=0, psym=4,symsize=Sym_Size $
  ,xrange=[x_min,x_max], xstyle=1  $
  ,yrange=[-25.,+0.], ystyle=1  $
  ;,xtitle='UT [hr]' $
  ,ytitle='MLat [deg]' $
, charsize = char_size, charthick = char_thick 

;3;mlt
plot, xplt0, zmlt,color=0, psym=4,symsize=Sym_Size $
  ,xrange=[x_min,x_max], xstyle=1  $
  ,yrange=[0.,+24.], ystyle=1  $
  ,xtitle='UT [hr]' $
  ,ytitle='MLT [hr]' $
, charsize = char_size, charthick = char_thick

output_png,flnm_png

save,/variables, filename='Ake_IPE'+cmpDate+'.'+version+'.sav'

print,'pro rd_akebono v1 finished'

end;pro rd_akebono v1
