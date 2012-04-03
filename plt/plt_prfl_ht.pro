;pro plot_profile_ht
pro plt_prfl_ht
fac_window=1.
sw_dif=0
sw_output2file=1

HOME_DIR='/home/Naomi.Maruyama/iper'
TEST='v55';43';
plot_UT =$
;16.00
24.002;v55
;23.07;v43
;23.83;v45
; freq_plot_hr = (60.*12.)/60. ;frequency of plotting in hr
freq_plot_hr=300./3600.;0.25000
;READ, freq_plot_hr,PROMPT='Enter frequency of plotting in hour:' 

title_hemi='NH'
;READ, title_hemi,PROMPT="Enter which hemisphere?: NH or SH" 

sw_fort=168L
;READ, sw_fort,PROMPT="Enter which fort?: 167 or 168"



sw_output2fil1e='PNG' ;NONE'
mp_plot=4-1L
lp_title=100 ;64;11 ;46
if ( lp_title eq 10 ) then $
  mlat_title='79.5' $
else if ( lp_title eq 11 ) then $
  mlat_title='78.4' $
else if ( lp_title eq 65 ) then $
  mlat_title='24.98' $
else if ( lp_title eq 64 ) then $
  mlat_title='25.26'  $
else if ( lp_title eq 100 ) then $
  mlat_title='16.664'  ;$
;if ( lp_title eq 41 ) then $
;  mlat_title='38'  $
;else if ( lp_title eq 42 ) then $
;  mlat_title='37'  $
;else if ( lp_title eq 43 ) then $
;  mlat_title='36'  $
;else if ( lp_title eq 44 ) then $
;  mlat_title='35'  $
;else if ( lp_title eq 45 ) then $
;  mlat_title='32'  $
;else if ( lp_title eq 46 ) then $


sw_debug=1L
plot_DIR=$
HOME_DIR+'/fig/'+TEST+'/'
;'../figures/discon/1dnewflipgrid/fort'+STRTRIM( string(sw_fort, FORMAT='(i3)'), 1)+'/'
n_file=1L

;high resolution
;if ( mlat_title eq '10' ) then  FLDIM =  265L  else $
;if ( mlat_title eq '31' ) then  FLDIM =  495L  else $
;if ( mlat_title eq '32' ) then  FLDIM = 1317L  else $
;if ( mlat_title eq '35' ) then  FLDIM = 1357L  else $
;if ( mlat_title eq '36' ) then  FLDIM = 1393L  else $  ;lp=43, mlat=35.8159
;if ( mlat_title eq '37' ) then  FLDIM = 1435L  else $  ;lp=42
;if ( mlat_title eq '38' ) then  FLDIM = 1473L  else $  ;lp=41
;if ( mlat_title eq '60' ) then  FLDIM = 2141L  else $
;if ( mlat_title ge '85' ) then  FLDIM = 4501L
FLDIM0=LONARR(n_file)
FLDIM0=[ $
; 161L ];lp=100 APEX grid v43
 401L ];lp=100 FLIP grid v55
;401L,401L,401L,401L]
;769L, 769L]
;769L, 695L, 651L, 843L];, 769L]
; 769L, 695L, 769L]
; 769L, 695L, 651L];, 769L]
; 769L, 695L,3831L, 3629L ]
;3831L, 3629L, 769L, 695L]
;3831L, 3629L, 3831L, 3629L]
;3831L, 3831L, 3831L, 3831L]  ;lp10, 3629 ;831L]
;3629L,3629L,3629L,3629L]        ; lp11
;431L, 431L, 431L, 431L]        ; lp65   
;433L, 433L, 433L, 433L]        ; lp64   
;433L, 431L,  433L, 431L ]
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
   open_fort168, n_file,LUN,sw_debug,title_hemi,lp_title,sw_fort,TEST

   while ( EOF(LUN[0]) eq 0 ) do begin

     plot_x[*,*,*,*] =-999999999.999999
     for i = 0, n_file-1  do begin
print,'i',i,'n_file',n_file
FLDIM=FLDIM0[i]
FLDIM_plot[i]=FLDIM

     Z0 = fltarr(FLDIM)
if ( sw_fort eq 168 ) then begin
   TNX0 = fltarr(FLDIM)
    UN0 = fltarr(FLDIM) ;cm/s
   NNO0 = fltarr(FLDIM)
 PHION0 = fltarr(FLDIM)
 VOP0   = fltarr(FLDIM) ; O+ field alinged velocity
   EHT0 = fltarr(FLDIM) ; EHT(3,J):e heating rate   
    TI0 = fltarr(3,FLDIM)
     N0 = fltarr(4,FLDIM)
SUMION0 = fltarr(3,FLDIM)  ;SUMION(1,7,J)+SUMION(2,4,J)+SUMION(2,5,J)
 xionn0  = fltarr(4,FLDIM)
 eqn2d0 = fltarr(FLDIM)
nplsrd0 = fltarr(FLDIM)

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

;print, 'i', i,FLDIM
if ( i eq 0 ) then sw_version = 1  else $
if ( i eq 1 ) then sw_version = 3  else $
if ( i eq 2 ) then sw_version = 3  else $
if ( i eq 3 ) then sw_version = 3

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
           read_fort168_3d, LUN[i],UT_hr0,LT_hr0,Z0,TNX0,UN0,NNO0,EHT0,TI0,N0,PHION0,VOP0,SUMION0,sw_version,mp_plot,mp_read,sw_debug   
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
     if ( sw_version eq 3 ) then $  ;yet 20110810another new version
       plot_x[i_window,  2,0:FLDIM-1,i] =   ALOG10 ( XIONN0[3-1,0:FLDIM-1] ) $ ;[he+]
     else $
       plot_x[i_window,  2,0:FLDIM-1,i] =   ALOG10 ( N0[4-1,0:FLDIM-1] ) ;[he+]
       plot_x[i_window,  3,0:FLDIM-1,i] =   ALOG10 ( N0[3-1,0:FLDIM-1] ) ;[Min+]

;        !.. Sum minor ions N+, NO+, O2+, N2+ for electron density at low altitudes
;         N(3,J)=XIONN(4,J)+XIONN(5,J)+XIONN(6,J)+XIONN(7,J)+XIONN(8,J)

;       plot_x[0,  4,0:FLDIM-1,i] =   ALOG10 ( NNO[   0:FLDIM-1] ) ;[NO]

i_window=1L
       plot_x[i_window,  0,0:FLDIM-1,i] =   TI0[1-1,0:FLDIM-1] ;Ti
       plot_x[i_window,  1,0:FLDIM-1,i] =   TI0[3-1,0:FLDIM-1] ;Te
;       plot_x[1,  2,0:FLDIM-1,i] =   TNX0[   0:FLDIM-1] ;Tn

;       plot_x[2,  0,0:FLDIM-1,i] =   VOP0[0:FLDIM-1] *1.E-01
;       plot_x[2,  1,0:FLDIM-1,i] =   UN0[ 0:FLDIM-1]
;       plot_x[2,  2,0:FLDIM-1,i] =   N0[1-1,0:FLDIM-1] * VOP0[0:FLDIM-1] * 1.E-06

       plot_x[2,  0,0:FLDIM-1,i] =   PHION0[0:FLDIM-1]
;       plot_x[3,  0,0:FLDIM-1,i] =   EHT0[  0:FLDIM-1]
       plot_x[3,  0,0:FLDIM-1,i] =   UN0[  0:FLDIM-1]*1.0E-2 ;cm-1-->m-s

print ,'check Un',MAX(UN0),MIN(UN0)
if ( sw_version ge 1 ) then begin
       plot_x[2,  1,0:FLDIM-1,i] =   SUMION0[0,0:FLDIM-1]  ;SUMION(1,7,J)
;       plot_x[3,  2,0:FLDIM-1,i] =   SUMION0[1,0:FLDIM-1] 
;       plot_x[3,  3,0:FLDIM-1,i] =   SUMION0[2,0:FLDIM-1]
endif 

;if ( sw_debug eq 1 ) then $

;if ( i ge 2 and i le 3 ) then $
;print, 'i=',i,'MAX phion',MAX(phion0[0:FLDIM-1]), MAX(eht0[0:FLDIM-1]) ;,MAX(sumion0[0,0:FLDIM-1]), MAX(sumion0[1,0:FLDIM-1]), MAX(sumion0[2,0:FLDIM-1])
endif else if ( sw_fort eq 167 ) then begin

       plot_x[0,  0,0:FLDIM-1,i] =   ALOG10( SL[   0:FLDIM-1] )

       plot_x[1,  0,0:FLDIM-1,i] =   GL[ 0:FLDIM-1] *180./!PI ;geolat
       plot_x[1,  1,0:FLDIM-1,i] =   SZA[0:FLDIM-1] *180./!PI ;sza

;print,'min',MIN(SZA[0:FLDIM-1] *180./!PI) ,' max',MAX(SZA[0:FLDIM-1] *180./!PI) 

       plot_x[2,  0,0:FLDIM-1,i] =   BM[0:FLDIM-1]             ;1,39E+00
       plot_x[2,  1,0:FLDIM-1,i] =   GR[ 0:FLDIM-1]*1.0E-3     ;-9.5E+2

       plot_x[3,  0,0:FLDIM-1,i] =   ALOG10( O0[  0:FLDIM-1] )
       plot_x[3,  1,0:FLDIM-1,i] =   ALOG10( H0[  0:FLDIM-1] + HE0[ 0:FLDIM-1]  )
       plot_x[3,  2,0:FLDIM-1,i] =   ALOG10( N20[ 0:FLDIM-1] + O20[ 0:FLDIM-1] )
       plot_x[3,  3,0:FLDIM-1,i] =   ALOG10( N4S0[ 0:FLDIM-1] )
   endif

;print, 'O+',MIN(plot_x[0,0,*,i]), MAX(plot_x[0,0,*,i])
;print, 'Ti',MIN(plot_x[1,0,*,i]), MAX(plot_x[1,0,*,i])

     endfor ;i = 0, n_file-1 do begin

;dbg20110804
;if ( UT_hr[0] ge 40.97 ) then  $
;print, ( (UT_hr[0]-UT_hr0_save) MOD freq_plot_hr ), (UT_hr[0]-UT_hr0_save),freq_plot_hr, UT_hr[0],UT_hr0_save  ;$
;else  GOTO, JUMP1



print,'UT_hr[0]',UT_hr[0],'plot_UT',plot_UT
     if ( UT_hr[0] ge plot_UT ) $
;AND  ( (UT_hr[0]-UT_hr0_save) MOD freq_plot_hr LT 0.0001 )$
 then begin


print,'plotting UT', ( (UT_hr[0]-UT_hr0_save) MOD freq_plot_hr ), (UT_hr[0]-UT_hr0_save), freq_plot_hr, UT_hr[0], UT_hr0_save

print,'plt_prfl',sw_output2file ;debug

;title_hemi='NH' ;requirement for j0/j1
;       profile_ht_3d $
       prfl_ht $
,plot_x,plot_y, title_hemi,mlat_title,ut_hr,lt_hr,plot_DIR,FLDIM_plot,mp_plot,sw_debug,sw_fort $
,sw_dif,sw_output2file,n_file,fac_window

BREAK ;exit from the while loop

;STOP
;JUMP1:
     endif ;( UT_hr[0] eq plot_UT ) then begin

 endwhile                       ; ( EOF(LUN168) ne 0 ) then begin

;JUMP1: 
for i = 0, n_file-1  do    FREE_LUN, LUN[i]
end ;pro plot_profile_ht
end ;pro plt_prfl_ht