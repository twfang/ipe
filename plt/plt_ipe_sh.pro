;20130828: copied from plt_ipe.pro to be called from shell script
;20131219: hinotori: Vartype in ctr_lon_lat needs to be zero to be able to output Ne. the output files are saved in ~/wamns/hinotori/. ncount at the end of the run should be input to read_hinotori.pro
;20130520UNDERCONSTRUCTION!!!include Tn for leslie
;tmp20121128 temporary o+ is assigned to plot_z(6)(n+) instead of 3 for faster debug molecular ions
;.20120305: renamed from plot_contour_3d_dif.pro-->plt_ipe.pro
;20120228: sw_dif: useed for comparison of two different runs
;20111205: sw_save=2 to save plotting time!!!
;include parallel plasma velocity to help the debug!!!
;20140813: execute from shell script
pro plt_ipe_sh
sw_debug=getenv('sw_debug')
sw_output2file=getenv('sw_output2file') ;1'PNG' ;0NONE';
sw_wam_ipe=getenv('sw_wam_ipe')
TEST=getenv('TEST')
TEST2='80';S';
;20131209: output to ascii file
sw_output2file_ascii=0;2L
f107=getenv('f107')
print,'f107=',f107
nday=getenv('nday')
print,'nday=',nday
sw_version_io = getenv('sw_version_io')
runDate = getenv('runDate')
print,'runDate=',runDate


   luntmp=100L
   luntmp1=101L
alt=getenv('alt')
print,' alt=', alt
if ( sw_output2file_ascii ge 1 ) then begin
  ; chr_title='F107='+STRTRIM( string(f107, FORMAT='(i3)'), 1)
  ; chr_title1=STRTRIM( string(alt, FORMAT='(F5.0)'), 1)+'km'
  ; chr_title2='NDAY='+STRTRIM( string(nday, FORMAT='(i3)'), 1)
;   flnmtmp='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/sed_nmf2'+chr_title1+'.'+chr_title2+'.'+chr_title+'.dat'
   flnmtmp='/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/sed_nmf2_20130317mp20lp35.dat'
   openw,luntmp,flnmtmp, /GET_LUN
   print, 'SED file created:',flnmtmp

endif ;( sw_output2file_ascii eq 1 ) then begin
;
n_read_max = getenv('n_read_max');
n_plt_max = getenv('n_plt_max')  ;for quick plot
;if n_plt_max eq 0 then n_plt_max = n_read_max - 1
;
n_plt_min = getenv('n_plt_min')  ;for quick plot
print,'n_read_max=',n_read_max,' n_plt_max=',n_plt_max
plot_UT    = getenv('plt_ut') 
plot_UT_end= getenv('plt_ut_end')
print,'plot_ut=', plot_ut,' plot_ut_end=', plot_ut_end 


sw_quickPlot = getenv('sw_quickPlot')
print, 'sw_quickplot ', sw_quickplot 
;20140117; plot every X hour
sw_hourly_plot = getenv('sw_hourly_plot')
print,' sw_hourly_plot', sw_hourly_plot
  plotXhr = getenv('plotXhr');0.2500
if sw_hourly_plot eq 1 then begin
  print,' plotXhr', plotXhr
endif

print,  'sw_grid=', getenv('sw_grid')
if ( getenv('sw_grid') eq 0 ) then $
  title_res='low' $
else if ( getenv('sw_grid') eq 1 ) then $
  title_res='low20120709' $
else if ( getenv('sw_grid') eq 2 ) then $
  title_res='2xdyn' $
else if ( getenv('sw_grid') eq 3 ) then $
  title_res='td20120709'



print,  'VarType=', getenv('VarType')
  VarType_min = FIX( getenv('VarType') )
  VarType_max = FIX( getenv('VarType') )
  VarType_step = 1L
if sw_debug eq 1 then  begin
  print, "check TYPE (1)",TYPENAME(vartype_min),' (2)',TYPENAME(vartype_max),' (3)',TYPENAME(vartype_step)
  print, "check vartype val",vartype_min,' (2)',vartype_max,' (3)varType_step=',vartype_step
endif





fac_window=getenv('fac_window')
;!!!CAUTION!!! plot_UT needs to be float (INT has limited digit!!!)

rpath = getenv('RPATH')
print,'rpath= ', rpath
rundir = getenv('rundir')
print,'rundir=',rundir
input_DIR0 = $
;rpath+'/'+rundir+'/'
;rpath+'/'
getenv('REFDIR')+'/'
print,'input_DIR0= ',input_DIR0

plot_type=FIX( getenv('plot_type') ) ;0L ;0:contour; 1:ht profile; 2:LT-LAT contour; 3:LON-LAT contour; 4:refilling: 5:psphere, 6:tec
print, 'plot_type=', plot_type
if plot_type eq 0 or plot_type eq 1 or plot_type eq 4 or plot_type eq 7 then begin
   mp_plot=FIX( getenv('mp_plot') ) ; longitude sector to plot
   print, 'mp_plot=', mp_plot
   mpstart=mp_plot
   mpstep=+1
   mpstop=mp_plot ;mpstart+mpstep
endif ;plot_type






;0:mag; 1:geo; 2:LT-maglat
sw_frame=FIX( getenv('sw_frame') )
print, 'sw_frame=',sw_frame
sw_dif=0L
sw_hr=0L
sw_3DJ=0L
sw_anim=0L
; freq_plot_hr = (60.*12.)/60. ;frequency of plotting in hr
output_freq=900
freq_plot_hr=output_freq/3600.
;READ, freq_plot_hr,PROMPT='Enter frequency of plotting in hour:'
sw_save=0L ;0:no save; 1:save; 2:restore
;sw_save option does not work for plot_type=1
;READ, sw_save,PROMPT='Enter switch for saving readin data 0or1:'
htstrt=400.
htstop=400.
htstep=200.
;;;!WARNING conflict with vartype_min/max for plot_type=0
Varstrt=0
Varstop=0

STOP_TIME='230406'
TEST0='trans'
title_test=TEST0+'.'+TEST  ;trans.'+TEST
title_hemi='eq';glb';SH' ;glb';eq';
version='3d'


fig_DIR=''
n_file=21L;
input_flnm=[$
 '','','','','','' $
,'','','','','','' $
,'','','','','','' $
,'','','']
input_DIR = input_flnm
input_DIR[*] = input_DIR0
input_DIR[1]='/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/grid/plt/';theia
LUN  = INTARR(n_file)
sw_LUN  = INTARR(n_file)
sw_lun[0]=getenv('sw_lun0') ;ut_rec
sw_lun[1]=1
sw_lun[2]=1 ;o+
sw_lun[6]=getenv('sw_lun6') ;h+
sw_lun[8]=getenv('sw_lun8') ;he+
sw_lun[9]=getenv('sw_lun9') ;n+
sw_lun[10]=getenv('sw_lun10') ;no+
sw_lun[11]=getenv('sw_lun11') ;o2+
sw_lun[12]=getenv('sw_lun12') ;n2+
sw_lun[13]=getenv('sw_lun13') ;o+(2D)
sw_lun[14]=getenv('sw_lun14') ;o+(2P)
sw_lun[3]=getenv('sw_lun3') ;Te
sw_lun[7]=getenv('sw_lun7') ;Ti
sw_lun[4]=getenv('sw_lun4')   ;vo+
sw_lun[5]=getenv('sw_lun5')   ;vexbup
sw_lun[15]=getenv('sw_lun15') ;vexbe
sw_lun[16]=getenv('sw_lun16') ;vexbth
sw_lun[17]=getenv('sw_lun2013') ;sunlon
;#nm20161128:plasma17 should be modified to plasma19 later...
sw_lun[18]=getenv('sw_lun19') ;sza
sw_lun[20]=getenv('sw_lun20') ;ipe_grid_plasma_params.dummy:george's new io
if ( sw_lun[20] eq 1 ) then  sw_lun[2]=0
if ( sw_dif eq 1 ) then begin
LUNq  = INTARR(n_file)
   input_flnmq =input_flnm
   input_DIRq =input_flnmq
   rundateq='20120217'
   TESTq='v46'
   title_testq=TEST0+'.'+TESTq
   input_DIRq[*]=$
TESTq+'/'
;rundateq+'.'+version+'.'+title_testq+'/';backup20120223/'
   input_DIRq[1]=input_DIR[1]
endif





plot_DIR=getenv('PPATH')+'/' 
print,' plot_DIR=', plot_DIR
if ( sw_frame eq 0 ) then $
title_frame='mag' $
else if ( sw_frame eq 1 ) then $
title_frame='geo'  $
else if ( sw_frame eq 2 ) then $
title_frame='lt'
filename_sav=plot_DIR+rundate+'_'+version+'.'+title_res+title_frame+'.sav'

if title_res eq 'low' then begin
  NLP=170L;low res
  NPTS2D=44438L ;low res
endif else if title_res eq 'low20120709' OR title_res eq 'td20120709' then begin
  NLP=170L;low res
  NPTS2D=44514L ;low res
endif else if title_res eq 'high' then  begin
  NLP=209L;high res
  NPTS2D=165717L ;high res
endif else if title_res eq '2xdyn' then  begin
  NLP=93L;
  NPTS2D=31287L ;high res
endif else if title_res eq 'dyn' then  begin
  NLP=45L;
  NPTS2D=15857L ;high res
endif

NMP=80L
ISPEC=9L
ISPEV=4L
MaxFluxTube=1115L ;=FLDIM
FLDIM=1115L;927L  ;???used when plot_type=1

  
  sunlons1 =  -1.92900
  UT_hr = 0.00D0
  UT_hr_save = fltarr(n_read_max)
;  LT_hr = fltarr(  NMP,NLP)
NPAR   = 6L;
;i should be aware of the memory limit!!!
if ( plot_type eq 0 ) or ( plot_type eq 2 ) or ( plot_type eq 4 )  or ( plot_type eq 7 ) then begin
;  plot_z = fltarr(n_read_max,NPAR, NMP,NPTS2D)
  plot_z = fltarr(n_read_max,NPAR,mpstop+1,NPTS2D) ;to save memory!!!
  if ( plot_type eq 0 )  or ( plot_type eq 7 )  then $
  plot_VEXB = fltarr(n_read_max,NMP,NLP)
endif

if ( sw_3DJ eq 1 ) then $
  je_3d=fltarr(2,NPTS2D,NMP)

if ( sw_hr eq 1 ) then $
  hrate=fltarr(7,NPTS2D,1) ;  hrate=fltarr(7,NPTS2D,NMP)


sw_read_wind = getenv('sw_read_wind')
;if ( sw_read_wind eq 1 ) then $
luntmp7=0L
luntmp3=0L
luntmp8=0L
luntmp9=0L
luntmp10=0L
luntmp11=0L
 Vn_ms1=fltarr(3,NPTS2D,NMP)
 on_m3 =fltarr(  NPTS2D,NMP)
 n2n_m3=fltarr(  NPTS2D,NMP)
 o2n_m3=fltarr(  NPTS2D,NMP)
 tn_k  =fltarr(  NPTS2D,NMP)
; Un_ms1=fltarr(MaxFluxTube,NLP,NMP)
XIONN_m3 =fltarr(ISPEC,NPTS2D,NMP)
XIONV_ms1=fltarr(1,NPTS2D,NMP) ;fltarr(ISPEV,NPTS2D,NMP)
TE_TI_k  =fltarr(3,NPTS2D,NMP)
VEXB=fltarr(NMP,NLP,3) ;before2011-11-16.v18 only meridional transport version
;after 2011-11-16.v18: VEXB=fltarr(2,NPTS2D,NMP) ;zonal transport version
;NHEAT_mks=dblarr(      NPTS2D,NMP)
;hrate_mks=dblarr(    6,NPTS2D,NMP)

if ( sw_dif eq 1 ) then begin
   XIONN_m3q =fltarr(3,NPTS2D,NMP) ;fltarr(ISPEC,NPTS2D,NMP)
   XIONV_ms1q=fltarr(1,NPTS2D,NMP) ;fltarr(ISPEV,NPTS2D,NMP)
   TE_TI_kq  =fltarr(3,NPTS2D,NMP)
   VEXBq=fltarr(NMP,NLP)
endif

n_read_max2=n_read_max
if ( plot_type eq 2 ) then begin
;  n_read_max2 = (133106-61106)/300 +1 ;=289
  plot_zz = fltarr(n_read_max2,NLP*2L)
  plot_xx = fltarr(n_read_max2,NLP*2L)
  plot_yy = fltarr(n_read_max2,NLP*2L)
endif

z_km =fltarr(NPTS2D) ;,NMP)
mlat_deg  =fltarr(NPTS2D) ;,NMP) ;GL_rad
JMIN_IN=lonarr(NLP)
JMAX_IS=lonarr(NLP)
glon_deg=fltarr(NPTS2D,NMP)
glat_deg=fltarr(NPTS2D,NMP)
sza_rad =fltarr(NPTS2D,NMP)
;glon=[ $
;  288.3167, 288.3167, 288.3167, 299.8955, 288.3167 $ ;0-4
;, 288.3167, 288.3167, 318.7912, 288.3167, 288.3167 $ ;5-9
;, 288.3167, 288.3167, 288.3167, 288.3167, 288.3167 $ ;10-14
;, 288.3167, 288.3167, 288.3167,   7.1758,  11.7843 $ ;15-19
;,  11.7843,  11.7843,  11.7843,  11.7843,  11.7843 $ ;20-24
;,  11.7843,  11.7843,  11.7843,  11.7843,  11.7843 $ ;25-29
;,  11.7843,  11.7843,  11.7843,  11.7843,  11.7843 $ ;30-34
;,  11.7843,  11.7843,  11.7843,  11.7843, 103.8842 $ ;35-39
;, 103.8842, 103.8842, 103.8842, 103.8842, 103.8842 $ ;40-44
;, 103.8842, 103.8842, 103.8842, 103.8842, 103.8842 $ ;45-49
;, 103.8842, 103.8842, 103.8842, 103.8842, 103.8842 $ ;50-54
;, 103.8842, 103.8842, 103.8842, 103.8842, 194.2077 $ ;55-59
;, 194.2077, 194.2077, 194.2077, 194.2077, 194.2077 $ ;60-64
;, 194.2077, 194.2077, 194.2077, 194.2077, 194.2077 $ ;65-69
;, 194.2077, 194.2077, 194.2077, 194.2077, 194.2077 $ ;70-74
;, 194.2077, 194.2077, 194.2077, 194.2077, 283.8234 ] ;75-79


if ( sw_save eq 2 ) then begin
  sw_save_old=sw_save
  plot_type_old=plot_type
sw_output2file_old=sw_output2file
  restore, filename=filename_sav
  print,'restoring the file=',filename_sav
  sw_save=sw_save_old
  plot_type=plot_type_old
sw_output2file=sw_output2file_old
endif



;need to debug when sw_save=2 
for ht_plot=340., 340., htstep  do begin
for VarType=1, 1  do begin

if ( sw_save eq 2 ) then begin
plot_zz[*,*]=0.0
plot_yy[*,*]=0.0
plot_xx[*,*]=0.0
endif


if ( plot_type eq 6 ) then begin
  NX=19L
  NY=121L
  TEC3d =fltarr(n_read_max,nx,nmp)
  glat3d=fltarr(n_read_max,nx,nmp)
  glon3d=fltarr(n_read_max,nx,nmp)
endif ;( plot_type eq 6 ) then begin

ncount=-1L
n_read=-1L
n_plt=-1L
;if ( sw_save lt 2 ) then $
;  check_value=EOF(LUN[1]) $
;else if ( sw_save eq 2 ) then $
;  check_value=n_read-n_read_max2

  ;while ( check_value le 0 ) do begin
  while ( n_read lt n_read_max-1 ) do begin
    n_read = n_read+1
    if ( sw_debug eq 1 ) then print,'n_read=',n_read,n_read_max2


;George's new restart io
    TimeStepMin=getenv('TimeStepMin')
    if sw_wam_ipe eq 0 then begin
       ;standalone IPE
       nSecMIN=0
       nMinMIN=0
       nSecTimeStamp = (n_read * TimestepMin * 60) + nSecMIN  ;[seconds]
       UT_hr = FLOAT(nSecTimeStamp)/3600.
       print,' nSecTimeStamp=', nSecTimeStamp,' n_read=',n_read,' nSecMin=',nSecMIN,' UT_hr=',UT_hr
       if nSecTimeStamp lt 10 then $
          TimeStamp4Plot =  'iter_0000000'+STRTRIM( string( nSecTimeStamp, FORMAT='(i1)'), 1) $
;dbg          TimeStamp4Plot =  'iter_000000'+STRTRIM( string( nSecTimeStamp, FORMAT='(i1)'), 1) $
       else if nSecTimeStamp lt 100 then $
          TimeStamp4Plot =  'iter_000000'+STRTRIM( string( nSecTimeStamp, FORMAT='(i2)'), 1)  $
       else if nSecTimeStamp lt 1000 then $
          TimeStamp4Plot =  'iter_00000'+STRTRIM( string( nSecTimeStamp, FORMAT='(i3)'), 1)  $
;dbg          TimeStamp4Plot =  'iter_00000'+STRTRIM( string( nSecTimeStamp, FORMAT='(i2)'), 1) 
       else if nSecTimeStamp lt 10000 then $
          TimeStamp4Plot =  'iter_0000'+STRTRIM( string( nSecTimeStamp, FORMAT='(i4)'), 1)   $
       else if nSecTimeStamp lt 100000 then $
          TimeStamp4Plot =  'iter_000'+STRTRIM( string( nSecTimeStamp, FORMAT='(i5)'), 1) 

    endif else if sw_wam_ipe eq 1 then begin
       ;WAM_IPE
       nHrMIN=0 ;###CHANGE!
       nHrTimeStamp = nHrMIN
;       runDate0=15
       runDate0=16
;       runDate0=17
;       runYear0=2013
       runYear0=2015
       runDate=STRTRIM( STRING( runYear0, FORMAT='(i4)'), 1)+'03'+STRTRIM( STRING( runDate0, FORMAT='(i2)'), 1)           
       nMinMIN=-TimeStepMin   ;0
       nMinTimeStamp = (n_read+1) * TimeStepMin + nMinMIN
       UT_hr = FLOAT(nMinTimeStamp)/60.
       print,' nMinTimeStamp=', nMinTimeStamp,' runDate=',runDate,' n_read=',n_read,' UT_hr=',UT_hr
       if nMinTimeStamp ge 60 then begin
         nHrTimeStamp = FIX( nMinTimeStamp/60 ) + nHrTimeStamp
         print,' nHrTimeStamp=', nHrTimeStamp
           if nHrTimeStamp ge 24 then begin
              nHrTimeStamp = ( FIX(nMinTimeStamp/60) MOD 24 )
;dbg20170825: temporary solution
              runDate0 = runDate0+1 
              runDate=STRTRIM( STRING( runYear0, FORMAT='(i4)'), 1)+'03'+STRTRIM( STRING( runDate0, FORMAT='(i2)'), 1)           
              print,' runDate=',runDate
           endif
           nMinTimeStamp = ( nMinTimeStamp MOD 60 )
        endif                   ;nMinTimeStamp
       print,' nMinTimeStamp=', nMinTimeStamp,' nHrTimeStamp=', nHrTimeStamp
       if nMinTimeStamp lt 10 then begin
          if nHrTimeStamp lt 10 then $
             TimeStamp4Plot = runDate+'0'+STRTRIM( STRING( nHrTimeStamp, FORMAT='(i2)'), 1)+'0'+STRTRIM( STRING( nMinTimeStamp, FORMAT='(i2)'), 1) $
          else if nHrTimeStamp lt 100 then $
             TimeStamp4Plot = runDate+STRTRIM( STRING( nHrTimeStamp, FORMAT='(i2)'), 1)+'0'+STRTRIM( STRING( nMinTimeStamp, FORMAT='(i2)'), 1)
          
       endif else if nMinTimeStamp lt 100 then begin
          if nHrTimeStamp lt 10 then $
             TimeStamp4Plot = runDate+'0'+STRTRIM( STRING( nHrTimeStamp, FORMAT='(i2)'), 1)+STRTRIM( STRING( nMinTimeStamp, FORMAT='(i2)'), 1) $
          else if nHrTimeStamp lt 100 then $
             TimeStamp4Plot = runDate+STRTRIM( STRING( nHrTimeStamp, FORMAT='(i2)'), 1)+STRTRIM( STRING( nMinTimeStamp, FORMAT='(i2)'), 1)
       endif                    ;nMinTimeStamp


    endif                       ;sw_wam_ipe
    print,' TimeStamp4Plot=', TimeStamp4Plot






;opening files
if ( sw_save le 1 ) then begin
   ;if ( n_read eq 0 ) OR ( sw_lun[20] eq 1 ) then begin
   if ( n_read eq 0 ) OR ( sw_version_io eq 1 ) then begin
      open_file,  input_DIR, LUN,version,input_flnm $
               ,  sw_3DJ,sw_hr, sw_lun,title_res $
               ,  sw_debug $
               , TimeStamp4Plot, n_read $
               , input_DIR0,luntmp7, luntmp3, luntmp8, luntmp9, luntmp10, luntmp11,sw_wam_ipe
      if n_read eq 0 then sw_lun[0:1]=0
   endif; n_read
endif ;if sw_save

;nm20170815 sw_dif does not work for the moment
if ( sw_dif eq 1 ) then $
   open_file,  input_DIRq, LUNq,version,input_flnmq $
               ,sw_3DJ,sw_hr, sw_lun,title_res $
               ,sw_debug





      if ( sw_save le 1 ) then begin

        if ( n_read eq 0 ) then  read_grid,LUN,JMIN_IN,JMAX_IS,Z_km,mlat_deg,sw_debug,glat_deg,glon_deg,title_res




        read_plasma_bin,LUN,UT_hr,XIONN_m3,XIONV_ms1,TE_TI_k,VEXB,sw_debug $
,sw_3DJ,je_3d,sw_hr,hrate, sw_dif, sw_lun $
,NMP,sunlons1 $
,sza_rad

  UT_hr_save[n_read]=UT_hr



  if ( sw_read_wind eq 1 ) then $
     read_wind, UT_hr $
                ,Vn_ms1,tn_k,on_m3,n2n_m3,o2n_m3 $
                , luntmp7, luntmp3, luntmp8, luntmp9, luntmp10, luntmp11 $
                , MaxFluxTube, nlp, nmp, JMIN_IN,JMAX_IS ,sw_debug, sw_version_io

if ( sw_debug eq 1 ) then  print,'after read_wind: check on_m3', MAX( on_m3 ), MIN( on_m3 )



  if ( sw_dif eq 1 ) then  begin
     read_plasma_bin,LUNq,UT_hrq,XIONN_m3q,XIONV_ms1q,TE_TI_kq,VEXBq,sw_debug $
,sw_3DJ,je_3d,sw_hr,hrate, sw_dif, sw_lun



if ( sw_dif eq 1 ) and ( UT_hr ne UT_hrq ) then begin
print,'!STOP INVALID UT! UThr', UT_hr,' UTq', UT_hrq
stop
endif

IF ( UT_hr lt plot_UT/3600. ) THEN CONTINUE


  endif
  
endif ;( sw_save eq 1 ) then begin

if ( plot_type eq 0 ) or ( plot_type eq 2 ) or ( plot_type eq 4 )  or ( plot_type eq 7 )  then begin
   ipts=0L
   if ( sw_save le 1 ) then begin


    for mp=mpstart,mpstop,mpstep do begin ;NMP-1 do begin

;0: Ne electron density[m-3]
      k=0L
      ;if ( sw_lun[2] eq 1 ) then begin
        for ipts=0L,NPTS2D-1L do begin 
          for jth=0,ISPEC-1L do begin
            plot_z[n_read,k,mp,ipts] = plot_z[n_read,k,mp,ipts] + XIONN_m3[jth,ipts,mp] 
          endfor;jth
        endfor ;ipts
      ;endif  ;sw_lun[2


;1 Te electron temperature
      jth=3-1
      k=1L
      if ( sw_lun[3] eq 1 ) then begin
        for ipts=0L,NPTS2D-1L do  plot_z[n_read,k,mp,ipts] = TE_TI_k[jth,ipts,mp]
      endif



;2 TO+ ion temperature
      jth=1-1
      k=2L                    
;t      if ( sw_lun[4] eq 1 ) then $
        for ipts=0L,NPTS2D-1L do  plot_z[n_read,k,mp,ipts] = XIONV_ms1[0,ipts,mp] ;dbg20140815 vo+
;       plot_z[n_read,k,mp,ipts] = TE_TI_k[jth,ipts,mp] 

;3-11  ion densities 
;for k=3,3+8 do begin
if ( sw_lun[2] eq 1 ) then begin
;for k=3,7 do begin
for k=3,5 do begin  ;dbg20140815
  jth=k-3
      for ipts=0L,NPTS2D-1L do $ 
       plot_z[n_read,k,mp,ipts] = XIONN_m3[jth,ipts,mp] 
endfor;k
endif ;( sw_lun[2] eq 1 ) then begin
;tmp20140130 temporary assign Vn_ms1 (field aligned wind) to k=5
; Vn_ms1=fltarr(3,NPTS2D,NMP)
if ( sw_read_wind eq 1 ) then begin
k=4L
jth=2-1
      for ipts=0L,NPTS2D-1L do $ 
;need to convert north-->southward as vo+ is positive southwrd
       plot_z[n_read,k,mp,ipts] = (-1.)*Vn_ms1[jth,ipts,mp] 
endif ;( sw_read_wind eq 1 ) then begin
;tmp20121128 temporary o+ is assigned to 6(n+) instead of 3 for faster debug molecular ions
;      for ipts=0L,NPTS2D-1L do $ 
;       plot_z[n_read,6,mp,ipts] = XIONN_m3[0,ipts,mp] 

;12:O+,velocity [m s-1]
k=12
;      for ipts=0L,NPTS2D-1L do $ 
;        plot_z[n_read,k,mp,ipts] = XIONV_ms1[1-1,ipts,mp]

;6:O+, flux [m-2 s-1]
;    for jth=0,0 do begin
;      plot_z[n_read,jth+6,mp,0:NPTS2D-1] = XIONN_m3[jth,0:NPTS2D-1,mp] * XIONV_ms1[jth,0:NPTS2D-1,mp]
;  endfor


;neutral heating rate eV/kg/s
if ( sw_hr eq 1 ) then begin 
  for jth=1,7 do begin
    plot_z[n_read,jth+1,mp,0:NPTS2D-1] = hrate[jth-1,0:NPTS2D-1,mp]
  endfor
endif
;    plot_z[5,   mp,0:NPTS2D-1] = NHEAT_mks[    0:NPTS2D-1,mp]
;    plot_z[6:11,mp,0:NPTS2D-1] = hrate_mks[0:5,0:NPTS2D-1,mp]

  if ( plot_type eq 0 ) or ( plot_type eq 7 ) then begin
        for lp=0,NLP-1 do begin
          midpoint = JMIN_IN(lp) + ( JMAX_IS(lp) - JMIN_IN(lp) )/2 -1
          plot_VEXB[n_read,mp,lp] = VEXB[mp,lp,2]           ;before 2011-11-16.v18:          
;after 2011-11-16.v18:          plot_VEXB[n_read,mp,lp] = VEXB[1,midpoint,mp]           
        endfor
  endif ;( plot_type eq 0 ) then begin

      endfor ;mp=0,NMP-1 do begin
  endif ;( sw_save eq 1 ) then begin

    if ( plot_type eq 2 ) then begin 
; LT-LAT plot
;dbg
;print,plot_type,'n_read',n_read
     ctr_LT_LAT $
  , mp_plot  $ 
  , JMIN_IN,JMAX_IS,Z_km,mlat_deg  $ 
  , plot_z   $
  , UT_hr_save, LT_hr, plot_DIR $
  , plot_zz,plot_xx,plot_yy,n_read $
  , VarType $
,ht_plot,sw_output2file $
,glon_deg,glat_deg,sw_frame $
, plot_UT
  endif ;plot_type eq 0


; lon-lat plot
    endif else if ( plot_type eq 3 ) then begin 

IF ( UT_hr lt plot_UT/3600. ) THEN CONTINUE
;if ( sw_debug eq 1 ) then  $
;print,'checkMOD=',(UT_hr MOD plotXhr),0.25,UT_hr,plotXhr
if  (sw_hourly_plot eq 1) AND ( (UT_hr MOD plotXhr) ge 0.010 ) then continue

n_plt = n_plt + 1
;dbg20170421 WAM IPE validation
;dbg20171114 temporary commented out to use the identical ctr routines!
;if sw_wam_ipe eq 1 then begin 
   ctr_lon_lat_quick_wam_ipe $
      , JMIN_IN,JMAX_IS,Z_km,mlat_deg  $ 
      , XIONN_m3, TE_TI_k $
      , XIONV_ms1 $
      , UT_hr, plot_DIR $
      , n_plt $
      , sw_output2file $
      ,glon_deg,glat_deg,sw_frame,fac_window, TEST $
      , sw_debug $
      , sw_output2file_ascii,luntmp,luntmp1,ncount $
      , Vn_ms1,tn_k,on_m3,n2n_m3,o2n_m3 $
      , n_plt_max,input_DIR0 $
      , alt,rundir $
      , VarType_max, VarType_min, VarType_step $
      , n_plt_min

;endif else if sw_wam_ipe eq 0 then begin 
;
;   if ( sw_quickplot eq 0 ) then $
;      ctr_lon_lat $
;      , JMIN_IN,JMAX_IS,Z_km,mlat_deg  $ 
;      , XIONN_m3, TE_TI_k $
;      , XIONV_ms1 $
;      , UT_hr, plot_DIR $
;      , n_read $
;      , sw_output2file $
;      ,glon_deg,glat_deg,sw_frame,fac_window, TEST $
;      , sw_debug $
;      , sw_output2file_ascii,luntmp,ncount $
;      , Vn_ms1,VEXB, sunlons1 $
;      , alt,rundir, LUN9001 $
;      , VarType_Min,VarType_Max, VarType_Step $
;   else if ( sw_quickplot eq 1 ) then $
;      ctr_lon_lat_quick $
;      , JMIN_IN,JMAX_IS,Z_km,mlat_deg  $ 
;      , XIONN_m3, TE_TI_k $
;      , XIONV_ms1 $
;      , UT_hr, plot_DIR $
;      , n_plt $
;      , sw_output2file $
;      ,glon_deg,glat_deg,sw_frame,fac_window, TEST $
;      , sw_debug $
;      , sw_output2file_ascii,luntmp,luntmp1,ncount $
;      , Vn_ms1,tn_k,on_m3,n2n_m3,o2n_m3 $
;      , n_plt_max,input_DIR0 $
;      , alt,rundir $
;      , VarType_max, VarType_min, VarType_step $
;      , n_plt_min;
;
;endif                           ;sw_wam_ipe eq 0 then begin 


;nm20161204 kitamura plot
    endif else if ( plot_type eq 5 ) then begin 

       plt_r_sza_ne $
  , JMIN_IN,JMAX_IS,Z_km, sza_rad $
  , xionn_m3, te_ti_k $
  ,sw_debug,mlat_deg,ut_hr,input_DIR0,plot_DIR

    endif else if ( plot_type eq 6 ) then begin 

;IF ( UT_hr lt plot_UT/3600. ) THEN CONTINUE


for mp=0,nmp-1 do begin
print,'plot_type=',plot_type,' mp=',mp

;interpolate Ne from 
ipdim=1115L
elden=fltarr(3,NPTS2D)
for i=1-1,NPTS2D-1 do begin
  for jth=0,ISPEC-1 do  elden[0,i]=elden[0,i]+XIONN_m3[jth,i,mp] 
  elden[1,i]=glat_deg[i,mp] 
  elden[2,i]=glon_deg[i,mp] 
endfor

;loadct,39
;xloadct
;stop
fxht_var=fltarr(3,nx,ny)
intrp2fixht $
,jmin_in,jmax_is,z_km $
,ut_hr $ 
,elden $
,ipdim,nlp,mlat_deg $
,fxht_var,fxht_x,fxht_y
;xloadct
;STOP
;calculate tec
TEC2d=fltarr(nx)
cal_tec $
,fxht_var,fxht_x,fxht_y $
,tec2d
;xloadct
tec3d[ n_read,0:nx-1,mp]=   tec2d[  0:nx-1]
print,'glat3d MIN=',MIN(fxht_var[1,0:nx-1,0]),'MAX',MAX(fxht_var[1,0:nx-1,0])
glat3d[n_read,0:nx-1,mp]=fxht_var[1,0:nx-1,0]
glon3d[n_read,0:nx-1,mp]=fxht_var[2,0:nx-1,0]
print,mp, MIN(tec2d[0:nx-1]), MAX(tec2d[0:nx-1])

;loadct,39
;lclr=(255/nmp) * mp
;print,'lclr',lclr
;if ( mp eq 0 ) then $
;  plot,  fxht_x,tec2d $
;,yrange=[ 12.,+15.], ystyle=1 $
;,xrange=[-60.,+60.], xstyle=1 $
;,color=lclr $
;else $ ;if ( mp eq 0 ) then $
;  oplot,  fxht_x,tec2d $
;,color=lclr ;$
;plt_tec, tec3d, fxht_x; , mlon???

endfor ;mp=0,nmp-1 do begin

    endif; plot_type eq 0

;dbg  endif ;( sw_save eq 1 ) then begin



;20130529: WHY these IF statements are not working!!!
;IF ( UT_hr lt plot_UT/3600. ) THEN CONTINUE
IF ( UT_hr gt plot_UT_end/3600. ) THEN BREAK
IF ( UT_hr ge plot_UT/3600. ) THEN BEGIN 

;
;



;
;if ( sw_debug eq 1 ) then  
;print, 'after plotting',UT_hr , plot_UT/3600. ;why this print sentence does not work???


;IF ( ( (UT_hr-UT_hr_save[0]) MOD freq_plot_hr ) LT 0.24 ) THEN BEGIN  ;20120223! need debug! does not work!!!


;glon_deg=fltarr(NPTS2D,NMP)
;jicamarca

if ( title_res eq 'low') OR ( title_res eq 'low20120709' ) OR ( title_res eq 'td20120709' )  then $
  lp=129L $;low
else if title_res eq '2xdyn' then $
  lp=70L $;2xdyn
else if title_res eq 'dyn' then $
  lp=34L ;dyn

midpoint= JMIN_IN[lp] + ( JMAX_IS[lp] - JMIN_IN[lp] ) /2

lt_hr2D=fltarr(NMP)
glon_deg2D=fltarr(NMP)
for mp=0,NMP-1 do begin
  lt_hr2D[mp]=UT_hr + glon_deg[midpoint,mp]/15.
  if ( lt_hr2D[mp] ge 24. ) then lt_hr2D[mp]=lt_hr2D[mp] MOD 24.
  glon_deg2D[mp]=glon_deg[midpoint,mp]
endfor

  if ( plot_type eq 0 ) then begin 
     if ( sw_debug eq 1 ) then     print, 'plotting contour: UT=',ut_hr




if( sw_debug eq 1 ) then  print,'plot_type',plot_type

if  (sw_hourly_plot eq 1) AND ( (UT_hr MOD plotXhr) ge 0.050 ) then continue

n_plt = n_plt + 1
if ( sw_quickplot eq 0 ) then $
     contour_plot_2d    $ 
  , JMIN_IN,JMAX_IS,Z_km,mlat_deg  $ 
  , plot_z,plot_VEXB,n_plt   $
  , UT_hr, plot_DIR, title_res $
,sw_debug, title_hemi,sw_anim,mpstart,mpstop,mpstep, lt_hr2D, fac_window $
  , sw_output2file, TEST $
,VarType_min $
,VarType_max $
,VarType_step $
;, input_DIR0,TEST,  TEST1, TEST2 $
else if ( sw_quickplot eq 1 ) then $
     ctr_lat_ht_quick    $ 
  , JMIN_IN,JMAX_IS,Z_km,mlat_deg  $ 
  , plot_z,plot_VEXB,n_read   $
  , UT_hr, plot_DIR, title_res $
,sw_debug, title_hemi,sw_anim,mpstart,mpstop,mpstep, lt_hr2D, fac_window $
  , sw_output2file $
,VarType_min $
,VarType_max $
,VarType_step $
,n_plt_max $
, input_DIR0,TEST,  TEST1, TEST2, glon_deg2D,rundir $
,n_plt_min ,tn_k,on_m3







  endif else if ( plot_type eq 1 ) then begin 
    print, 'plotting ht profile: UT=',ut_hr

plot_type_prof=0L ;0:densities; 1:temperatures
title_hemi='NH'
n_file_plot =1L
ut_hr_plot = fltarr(n_file_plot)
ut_hr_plot[0:n_file_plot-1] = UT_hr
lt_hr_plot = fltarr(n_file_plot)
 FLDIM_max = FLDIM  ;for mlat_title='38' lp=41
if ( sw_debug eq 1 ) then  print, 'FLDIM_max', FLDIM_max,'n_file_plot',n_file_plot
 plot_y = fltarr(  FLDIM_max,n_file_plot)
plot_type_max=4L
k_species = 4L
plot_x = fltarr(plot_type_max, k_species,FLDIM_max,n_file_plot)


;dbg
;d print,'n_file_plot',n_file_plot
;d print,'size(plot_x)',size(plot_x)
;d print,'size(plot_y)',size(plot_y)

     plot_x[*,*,*,*] =-999999999.999999
 FLDIM_plot=LONARR(n_file_plot)
for   mp_plot0=mp_plot,mp_plot,1  do begin
  lp_plot0= 1-1L

for i_file = 0,n_file_plot-1 do begin

if (i_file eq 0 ) then begin
  lp_plot=lp_plot0
;d  mp_plot=mp_plot0
endif else if (i_file eq 1 ) then begin
  lp_plot=136-1 ;lp_plot0-1L
;d  mp_plot=mp_plot0
endif else if (i_file eq 2 ) then begin
  lp_plot=lp_plot0+0L
;d  mp_plot=mp_plot0
endif else if (i_file eq 3 ) then begin
  lp_plot=lp_plot0+1L
;d  mp_plot=mp_plot0
endif else if (i_file eq 4 ) then begin
  lp_plot=lp_plot0+2L
;d  mp_plot=mp_plot0
endif

;  if (i_file eq 0 ) then $
;    lp_plot=10-1L $  ;idl conversion
;  else if (i_file eq 1 ) then $
;    lp_plot=11-1L $   ;idl conversion
;  else if (i_file eq 2 ) then $
;    lp_plot=12-1L $   ;idl conversion
;  else if (i_file eq 3 ) then $
;    lp_plot=9-1L   

print,'i_file',i_file,'lp_plot',lp_plot,'mp_plot0',mp_plot0

  lp_title=lp_plot+1 ;=42
;  if ( lp_title eq 40 ) then $
;    mlat_title='39'  $  ;???
;  else if ( lp_title eq 41 ) then $
;    mlat_title='38'  $
;  else if ( lp_title eq 42 ) then $
;    mlat_title='37'  $
;  else if ( lp_title eq 43 ) then $
;    mlat_title='36'  $
;  else if ( lp_title eq 44 ) then $
;    mlat_title='35'  $
;  else if ( lp_title eq 45 ) then $
;    mlat_title='32'  $
;  else if ( lp_title eq 46 ) then $
;    mlat_title='31'  $
;  else if ( lp_title eq 47 ) then $
;    mlat_title='29'   



  in=JMIN_IN[LP_plot]-1  ;idl conversion
  is=JMAX_IS[LP_plot]-1  ;idl conversion
  FLDIM_plot[i_file] = is-in+1L
  midpoint = JMIN_IN(lp_plot) + ( JMAX_IS(lp_plot) - JMIN_IN(lp_plot) )/2 -1
;  lt_hr_plot[i_file] = LT_hr ;[mp_plot,lp_plot]
  lt_hr_plot[i_file] = ut_hr + glon_deg[midpoint,mp_plot]/15.

print,'FLDIM=',FLDIM_plot[i_file],in,is,' mlat[deg]',mlat_deg[in]
mlat_title=STRTRIM( string(mlat_deg[in], FORMAT='(f8.3)'), 1)

factor=1.0E-11
  if( sw_dif eq 0 ) then begin  
     i_window=0L
     plot_x[i_window,0,0:FLDIM_plot[i_file]-1,i_file] =$
;   ALOG10 ( XIONN_m3[1-1,in:is,MP_plot] * 1.0E-6 ) ;[o+]m-3 --> cm-3
XIONN_m3[1-1,in:is,MP_plot] * factor ;[o+]m-3 X10^11
     plot_x[i_window,1,0:FLDIM_plot[i_file]-1,i_file] =$
;   ALOG10 ( XIONN_m3[2-1,in:is,MP_plot] * 1.0E-6 ) ;[h+]
XIONN_m3[2-1,in:is,MP_plot] * factor ;[h+]
     plot_x[i_window,2,0:FLDIM_plot[i_file]-1,i_file] =$
;   ALOG10 ( XIONN_m3[3-1,in:is,MP_plot] * 1.0E-6 ) ;[he+]
XIONN_m3[3-1,in:is,MP_plot] * factor ;[he+]
     plot_x[i_window,3,0:FLDIM_plot[i_file]-1,i_file] =$
;   ALOG10 ( XIONN_m3[4-1,in:is,MP_plot] * 1.0E-6 ) ;[N+]
XIONN_m3[4-1,in:is,MP_plot] * factor ;[N+]
     plot_x[i_window,0,0:FLDIM_plot[i_file]-1,i_file] =   ALOG10 ( XIONN_m3[1-1,in:is,MP_plot] * 1.0E-6 ) ;STORM [o+]m-3 --> cm-3
     plot_x[i_window,1,0:FLDIM_plot[i_file]-1,i_file] =   ALOG10 ( XIONN_m3[1-1,in:is,MP_plot] * 1.0E-6 ) ;Q [o+]
     plot_x[i_window,2,0:FLDIM_plot[i_file]-1,i_file] =   ALOG10 ( XIONN_m3[2-1,in:is,MP_plot] * 1.0E-6 ) ;S [h+]
     plot_x[i_window,3,0:FLDIM_plot[i_file]-1,i_file] =   ALOG10 ( XIONN_m3[2-1,in:is,MP_plot] * 1.0E-6 ) ;Q [h+]
     i_window=1L
     plot_x[i_window,0,0:FLDIM_plot[i_file]-1,i_file] =   TE_TI_k[1-1,in:is,MP_plot] ;STORM Ti
     plot_x[i_window,1,0:FLDIM_plot[i_file]-1,i_file] =   TE_TI_k[1-1,in:is,MP_plot] ;QUIET Ti
     plot_x[i_window,2,0:FLDIM_plot[i_file]-1,i_file] =   TE_TI_k[3-1,in:is,MP_plot] ;STORM Te
     plot_x[i_window,3,0:FLDIM_plot[i_file]-1,i_file] =   TE_TI_k[3-1,in:is,MP_plot] ;QUIET Te
  endif 

  plot_y[    0:FLDIM_plot[i_file]-1,i_file] =                  Z_km[  in:is]
;endif
endfor                     ;i_file = 0,n_file_plot-1 do begin

sw_fort=168L
;20120305:     profile_ht_3d $
     prfl_ht $
        ,plot_x,plot_y, title_hemi,mlat_title,ut_hr_plot,lt_hr_plot $
        ;,plot_type_prof
        ,plot_DIR,FLDIM_plot,mp_plot0,sw_debug,sw_fort $
        ,sw_dif,sw_output2file,n_file_plot,fac_window,TEST,rundir
endfor   ;mp_plot0=1-1,nmp-1,10  do begin     


;nm20170208 plot plasmasphere2D for akebono comparison
  endif else if ( plot_type eq 7 ) then begin 
print,n_read,ut_hr, sw_anim,mp_plot
       ctr_plt_psphere    $ 
          , JMIN_IN,JMAX_IS,Z_km,mlat_deg  $ 
          , plot_z,plot_VEXB,n_read   $
          , UT_hr, plot_DIR, title_res,rundate,title_test,sw_debug, title_hemi,sw_anim,mp_plot,  fac_window $
          , sw_output2file

  endif ;else if ( plot_type eq 1 ) then begin 
;BREAK ;exit from while loop

  if ( UT_hr ge plot_UT_end/3600. ) and ( plot_type eq 4 ) then begin

;refilling
    plt_refil, mlat_deg, JMIN_IN,JMAX_IS, plot_z,mp_plot,n_read_max,ut_hr_save,fac_window, plot_UT,plot_UT_end, TEST,plot_DIR,sw_debug

  endif ;( plot_type eq 4 ) then 

endif ;( UT_hr eq plot_UT ) then begin
endwhile                            ; ( EOF(LUN[0]) ne 0 ) then begin
endfor ;VarType=1,7 do begin
endfor ;ht_plot=200.00,400.00, 100.00  do begin

if ( sw_save eq 1 ) then  begin
  print,'saving to a file=',filename_sav
  save, filename=filename_sav
  print,'saving finished'
endif

if ( sw_save le 1 ) then begin 
  for i = 0, n_file-1  do   if ( sw_lun[i] eq 1 ) then   FREE_LUN, LUN[i]
  if ( sw_dif eq 1 ) then $
    for i = 0, n_file-1  do   if ( sw_lun[i] eq 1 ) then   FREE_LUN, LUNq[i]
endif

if ( sw_output2file_ascii eq 1 ) then begin
  FREE_LUN, luntmp ;openw,luntmp,flnmtmp, /GET_LUN
  print, 'ncount', ncount
endif

print,'plt_ipe_sh: finished successfully!'
end ;pro plt_ipe
