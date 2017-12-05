;20140610: sw_dif: underconstruction
;pro plot_profile_ht
pro plt_prfl_ht
todaysDate='20171129'
n_read_max=302L
fac_window=1.5
sw_dif=0
sw_output2file=1L;1L
VTEST='2015031500'
TEST0=$
;'20160727fricHeat1'
;'20160916eldyn'
;'20171116IPEOptimization'
;'20171109raw_high_latWamIpe'+VTEST
'20171126raw_high_latWamIpe'+VTEST
HOME_DIR=$
;'/scratch3/NCEPDEV/swpc/save/Tzu-Wei.Fang/IPE/trunk_electrodynamics_44514_feedback/run/'
;'/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/mergeEldyn20170127Stefan/ipe/run/'
;'/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/stmp'+VTEST+'/Naomi.Maruyama/'
'/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/stmp'+VTEST+'/Naomi.Maruyama/'
;'/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/20171025IpeOptimization/ipe/run/'
;CHANGE!!!
rundir=$
'variablef107kp'
;'new_compsets'
;'1510852798_ipe_theia_intel_parallel_464'
plot_UT =0.0;529200./3600.;2.0
;frequency of plotting in hr
freq_plot_hr=1.00000;0.01666666 ;60./3600.
;READ, freq_plot_hr,PROMPT='Enter frequency of plotting in hour:' 

title_hemi='SH';
;READ, title_hemi,PROMPT="Enter which hemisphere?: NH or SH" 

sw_fort=167L;168L
;READ, sw_fort,PROMPT="Enter which fort?: 167 or 168"


mp_plot=5-1L ;<--read from the 167 file!

n_file=1L
FLDIM0=LONARR(n_file)

;CHANGE!!!
;if ( rundir eq 'ipe_80_17352dbg' ) then begin
;  lp_title    =44-1L ;<--read from the 167 file!
;   mlat_title ='34.37' ;<--read from 167 file! =GL/!PI*180. lp=44
;   FLDIM0     =[ 395L ] ;<--read from the 167 file! lp=44
;endif else if ( rundir eq 'ipe_80_7459dbg' ) then begin
;  lp_title    =48-1L
;  mlat_title  ='29.67' ;<--read from 167 file! =GL/!PI*180. lp=48
;  FLDIM0     =[ 313L ] ;<--read from the 167 file! lp=48
;endif
;mp=10-1, lp=14-1
lp_title    =103-1L
mlat_title  ='-14.0';5.84';82.4001';76.2';66.13496' ;<--read from 167 file! =GL/!PI*180. lp=48
FLDIM0     =[ 157L ] ;<--read from the 167 file! lp=48

sw_debug=0L
plot_DIR=$
'/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/fig/Naomi.Maruyama/'+todaysDate+'/prfl/test'+VTEST+'/'+title_hemi+'/'
;'/scratch3/NCEPDEV/stmp2/Naomi.Maruyama/eldyn/fig/'+rundir+'/prfl/'


  LUN  = INTARR(n_file)
 UT_hr = fltarr(n_file)
 LT_hr = fltarr(n_file)
FLDIM_max=MAX(FLDIM0)
    plot_y = fltarr(      FLDIM_max,n_file)
plot_type_max=4L
k_species = 4L
    plot_x = fltarr(plot_type_max, k_species, FLDIM_max,n_file)
FLDIM_plot=LONARR(n_file)

; temporary array just for reading
 UT_hr0 = 0.0
 UT_hr0_save = -99999.
 LT_hr0 = 0.0


;FLDIMphil=401
;opening files
   open_fort168, n_file,LUN,sw_debug,title_hemi,lp_title,sw_fort,rundir,TEST0,HOME_DIR


   n_read=-1L
   while ( EOF(LUN[0]) eq 0 ) do begin
      n_read=n_read+1

      if n_read gt n_read_max then BREAK ;exit from while EOF(LUN0) loop


      plot_x[*,*,*,*] =-999999999.999999
      for i = 0, n_file-1  do begin
         print,'i=',i,'n_file=',n_file
         FLDIM=FLDIM0[i]

         if ( n_read eq 0 ) then  print,'FLDIM=',FLDIM
         FLDIM_plot[i]=FLDIM
         
         Z0 = fltarr(FLDIM)
         if ( sw_fort eq 168 ) then begin
            TNX0 = fltarr(FLDIM)
            UN0 = fltarr(FLDIM) ;cm/s
            NNO0 = fltarr(FLDIM)
            PHION0 = fltarr(FLDIM)
            VOP0   = fltarr(FLDIM) ; O+ field alinged velocity
            EHT0 = fltarr(FLDIM)   ; EHT(3,J):e heating rate   
            TI0 = fltarr(3,FLDIM)
            N0 = fltarr(4,FLDIM)
            SUMION0 = fltarr(3,FLDIM) ;SUMION(1,7,J)+SUMION(2,4,J)+SUMION(2,5,J)
            xionn0  = fltarr(4,FLDIM)
            eqn2d0 = fltarr(FLDIM)
            nplsrd0 = fltarr(FLDIM)
            
            if ( sw_dif eq 1 ) then begin
               N1 = fltarr(4,FLDIM)
               TI1 = fltarr(3,FLDIM)
               PHION1 = fltarr(FLDIM)
               UN1 = fltarr(FLDIM) ;cm/s
               SUMION1 = fltarr(3,FLDIM) ;SUMION(1,7,J)+SUMION(2,4,J)+SUMION(2,5,J)
            endif

         endif else if ( sw_fort eq 167 ) then begin
            SL = fltarr(FLDIM)
            GL = fltarr(FLDIM)
            BM = fltarr(FLDIM)
            GR = fltarr(FLDIM)
            SZA = fltarr(FLDIM)
            O0 = fltarr(FLDIM)
            H0 = fltarr(FLDIM)
            N20 = fltarr(FLDIM)
            O20 = fltarr(FLDIM)
            HE0 = fltarr(FLDIM)
            N4S0 = fltarr(FLDIM)
         endif


         if ( i eq 0 ) then sw_version = 1  $
         else if ( i eq 1 ) then sw_version = 3  $
         else if ( i eq 2 ) then sw_version = 3  $
         else if ( i eq 3 ) then sw_version = 3


;       if ( i ge 0 and i le 1) then $
;if ( sw_fort eq 168L ) then $
;         read_fort168, LUN[i],UT_hr0,LT_hr0,Z0,TNX0,UN0,NNO0,EHT0,TI0,N0,PHION0,VOP0,SUMION0,sw_version,sw_debug $
;         ,xionn0, eqn2d0, nplsrd0 $
         if ( sw_fort eq 167L ) then $
            read_fort167, LUN[i],UT_hr0,LT_hr0,Z0,SL,GL,BM,GR,SZA,O0,H0,N20,O20,HE0,N4S0,sw_debug,title_hemi $

;       else if ( i ge 1 ) then begin 

;         if ( i eq 1 ) then begin 
;           mp_read=4-1L
;           read_fort168_3d, LUN[i],UT_hr0,LT_hr0,Z0,TNX0,UN0,NNO0,EHT0,TI0,N0,PHION0,VOP0,mp_plot,mp_read,sw_debug   
;       else if ( i ge 2 and i le 3) then begin 
         else if ( sw_fort eq 168L ) then begin
            mp_read=1-1L
            read_fort168_3d, LUN[i],UT_hr0,LT_hr0,Z0,TNX0,UN0,NNO0,EHT0,TI0,N0,PHION0,VOP0,SUMION0,sw_version,mp_plot,mp_read,sw_debug,n_read   
            if ( sw_dif eq 1 ) then  read_fort168_3d, LUN[i],UT_hr0,LT_hr0,Z0,TNX0,UN1,NNO0,EHT0,TI1,N1,PHION1,VOP0,SUMION1,sw_version,mp_plot,mp_read,sw_debug    
endif
;         endif ;else if ( i eq 2 ) then begin 
;       endif ;else if ( i ge 1 ) then begin 

;this is necessary if reading Phil's output
;else $
;if ( i eq 1 ) then $
;read_fort168phil, LUN[i],UT_hr0,LT_hr0,Z0,TNX0,UN0,NNO0,EHT0,TI0,N0,PHION0,sw_debug

         UT_hr[i] = UT_hr0
         if ( i eq 0 and UT_hr0_save le -99999. ) then  begin
            UT_hr0_save = UT_hr0
            print, 'UT_hr0_save=',UT_hr0_save , UT_hr0
         endif

         LT_hr[i] = LT_hr0
         plot_y[    0:FLDIM-1,i] =    Z0[    0:FLDIM-1]
         
         if ( sw_fort eq 168 ) then begin
            i_window=0L
            plot_x[i_window,  0,0:FLDIM-1,i] =   ALOG10 ( N0[1-1,0:FLDIM-1] ) ;[o+]
            plot_x[i_window,  1,0:FLDIM-1,i] =   ALOG10 ( N0[2-1,0:FLDIM-1] ) ;[h+]
            if ( sw_version eq 3 ) then $ ;yet 20110810another new version
               plot_x[i_window,  2,0:FLDIM-1,i] =   ALOG10 ( XIONN0[3-1,0:FLDIM-1] ) $ ;[he+]
            else $
               plot_x[i_window,  2,0:FLDIM-1,i] =   ALOG10 ( N0[4-1,0:FLDIM-1] ) ;[he+]
            plot_x[i_window,  3,0:FLDIM-1,i] =   ALOG10 ( N0[3-1,0:FLDIM-1] )    ;[Min+]
            
;        !.. Sum minor ions N+, NO+, O2+, N2+ for electron density at low altitudes
;         N(3,J)=XIONN(4,J)+XIONN(5,J)+XIONN(6,J)+XIONN(7,J)+XIONN(8,J)

;       plot_x[0,  4,0:FLDIM-1,i] =   ALOG10 ( NNO[   0:FLDIM-1] ) ;[NO]
            if ( sw_dif eq 1 ) then begin 
               plot_x[i_window,  0,0:FLDIM-1,i] =   ALOG10 ( N1[1-1,0:FLDIM-1] ) - ALOG10 ( N0[1-1,0:FLDIM-1] ) ;[o+]
               plot_x[i_window,  1,0:FLDIM-1,i] =   ALOG10 ( N1[2-1,0:FLDIM-1] ) - ALOG10 ( N0[2-1,0:FLDIM-1] ) ;[h+]
;     if ( sw_version eq 3 ) then $  ;yet 20110810another new version
;       plot_x[i_window,  2,0:FLDIM-1,i] =   ALOG10 ( XIONN0[3-1,0:FLDIM-1] )  $ ;
;     else $
               plot_x[i_window,  2,0:FLDIM-1,i] =   ALOG10 ( N1[4-1,0:FLDIM-1] ) - ALOG10 ( N0[4-1,0:FLDIM-1] ) ;[he+]
               plot_x[i_window,  3,0:FLDIM-1,i] =   ALOG10 ( N1[3-1,0:FLDIM-1] ) - ALOG10 ( N0[3-1,0:FLDIM-1] ) ;[Min+]
            endif               ;( sw_dif eq 1 ) then 

            i_window=1L
            plot_x[i_window,  0,0:FLDIM-1,i] = TI0[1-1,0:FLDIM-1] ;Ti
            plot_x[i_window,  1,0:FLDIM-1,i] = TI0[3-1,0:FLDIM-1] ;Te
            plot_x[i_window,  2,0:FLDIM-1,i] = TNX0[   0:FLDIM-1] ;Tn

            if ( sw_dif eq 1 ) then begin
               plot_x[i_window,  0,0:FLDIM-1,i] =   TI1[1-1,0:FLDIM-1] - TI0[1-1,0:FLDIM-1] ;Ti
               plot_x[i_window,  1,0:FLDIM-1,i] =   TI1[3-1,0:FLDIM-1] - TI0[3-1,0:FLDIM-1] ;Te
            endif


            i_window=2L
            plot_x[i_window,  0,0:FLDIM-1,i] =   ALOG10( PHION0[0:FLDIM-1] )

            if ( sw_dif eq 1 ) then $
               plot_x[i_window,  0,0:FLDIM-1,i] =   PHION1[0:FLDIM-1] - PHION0[0:FLDIM-1]
;       plot_x[2,  0,0:FLDIM-1,i] =   VOP0[0:FLDIM-1] *1.E-01
;       plot_x[2,  2,0:FLDIM-1,i] =   N0[1-1,0:FLDIM-1] * VOP0[0:FLDIM-1] * 1.E-06

            if ( sw_version ge 1 ) then begin
               plot_x[i_window,  1,0:FLDIM-1,i] = ALOG10(  SUMION0[0,0:FLDIM-1] );SUMION(1,7,J)

               if ( sw_dif eq 1 ) then $
                  plot_x[i_window,  1,0:FLDIM-1,i] =   SUMION1[0,0:FLDIM-1] - SUMION0[0,0:FLDIM-1] ;SUMION(1,7,J)
            endif 
            


            i_window=3L            
            plot_x[i_window,  0,0:FLDIM-1,i] =   UN0[  0:FLDIM-1]*1.0E-2 ;cm-1-->m-s

;       plot_x[3,  2,0:FLDIM-1,i] =   SUMION0[1,0:FLDIM-1] 
;       plot_x[3,  3,0:FLDIM-1,i] =   SUMION0[2,0:FLDIM-1]
            
            if ( sw_dif eq 1 ) then $
               plot_x[3,  0,0:FLDIM-1,i] =   UN1[  0:FLDIM-1]*1.0E-2 - UN0[  0:FLDIM-1]*1.0E-2 ;cm-1-->m-s



;if ( sw_debug eq 1 ) then $


         endif else if ( sw_fort eq 167 ) then begin

            plot_x[0,  0,0:FLDIM-1,i] =   ALOG10( SL[   0:FLDIM-1] )
            
            plot_x[1,  0,0:FLDIM-1,i] =   GL[ 0:FLDIM-1] *180./!PI ;geolat
            plot_x[1,  1,0:FLDIM-1,i] =   SZA[0:FLDIM-1] *180./!PI ;sza
            


            plot_x[2,  0,0:FLDIM-1,i] =   BM[0:FLDIM-1]    ;1,39E+00
            plot_x[2,  1,0:FLDIM-1,i] =   GR[ 0:FLDIM-1]*1.0E-3 ;-9.5E+2
            
            plot_x[3,  0,0:FLDIM-1,i] =   ALOG10( O0[  0:FLDIM-1] )
            plot_x[3,  1,0:FLDIM-1,i] =   ALOG10( H0[  0:FLDIM-1] + HE0[ 0:FLDIM-1]  )
            plot_x[3,  2,0:FLDIM-1,i] =   ALOG10( N20[ 0:FLDIM-1] + O20[ 0:FLDIM-1] )
            plot_x[3,  3,0:FLDIM-1,i] =   ALOG10( N4S0[ 0:FLDIM-1] )
         endif



      endfor                    ;i = 0, n_file-1 do begin




      if ( sw_debug eq 1 ) then  $
print,'before plotting: MOD=', ( (UT_hr[0]-UT_hr0_save) MOD freq_plot_hr ),' diffut=', (UT_hr[0]-UT_hr0_save), ' freq=',freq_plot_hr, ' UT=',UT_hr[0], ' UT0=', UT_hr0_save


      IF (  ((UT_hr[0]-UT_hr0_save) MOD freq_plot_hr) GE 0.002 ) then begin
         if ( sw_debug eq 1 ) then $
 print,'skip plotting'
         continue
      ENDIF

      if ( UT_hr[0] ge plot_UT )  then begin


         if ( sw_debug eq 1 ) then  $
print,'start plotting: MOD=', ( (UT_hr[0]-UT_hr0_save) MOD freq_plot_hr ),' diffut=', (UT_hr[0]-UT_hr0_save),' freq=', freq_plot_hr, ' UT=',UT_hr[0],' UT0=', UT_hr0_save




;       profile_ht_3d $
       prfl_ht $
,plot_x,plot_y, title_hemi,mlat_title,ut_hr,lt_hr,plot_DIR,FLDIM_plot,mp_plot,sw_debug,sw_fort $
,sw_dif,sw_output2file,n_file,fac_window, TEST0,rundir



;JUMP1:
     endif ;( UT_hr[0] eq plot_UT ) then begin

 endwhile                       ; ( EOF(LUN168) ne 0 ) then begin
print,'endwhile: n_read=',n_read,'n_read_max=',n_read_max

;JUMP1: 
for i = 0, n_file-1  do    FREE_LUN, LUN[i]
print, 'plt_prfl_ht finished successfully!'
end ;pro plot_profile_ht

