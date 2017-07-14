;20121025: make sure to prepare the grid file (binary) in
;ipe/plt/grd/plasma_grid.XXX
pro read_grid,LUN,JMIN_IN,JMAX_IS,Z_km,mlat_deg,sw_debug,glat_deg,glon_deg,title_res

plot_type=getenv('plot_type')

;20121025: .save did not save too much time as compared to the
;original binary file.....
sw_save_grid=0
filename_grid_sav='grd/plasma_grid.'+title_res+'.sav'

if ( sw_save_grid eq 2 ) then  begin
  LUN_old=LUN
  sw_debug_old=sw_debug
  title_res_old=title_res

  restore, filename=filename_grid_sav
  print,'restoring the grid file=',filename_grid_sav

  LUN=LUN_old
  sw_debug=sw_debug_old
  title_res=title_res_old
  RETURN
endif

     mp=3L

if sw_debug eq 1 then  print,'size jmin',SIZE(jmin_in)
if sw_debug eq 1 then  print,'size jmax',SIZE(jmax_is)




if sw_debug eq 1 then  print,'size z_km',SIZE(Z_km)


if sw_debug eq 1 then  print,'size mlat_deg ',SIZE(mlat_deg)

NMP0=0L
NMP1=0L
NMP_all=0L
NLP_all=0L
NPTS2D_dum=0L
; read plasma_grid
if ( title_res eq 'dyn' ) OR ( title_res eq 'low20120709' ) OR ( title_res eq '2xdyn' ) $
;OR ( title_res eq 'td20120709' ) $
then begin
   readu, LUN[1], NMP0
if sw_debug eq 1 then    print,'NMP0',NMP0
   readu, LUN[1], NMP1
if sw_debug eq 1 then    print,'NMP1',NMP1
endif

    readu, LUN[1], NMP_all
if sw_debug eq 1 then print,'NMP_all',NMP_all
    readu, LUN[1], NLP_all
if sw_debug eq 1 then print,'NLP_all',NLP_all
    readu, LUN[1], NPTS2D_dum
if sw_debug eq 1 then print,'NPTS2D_dum',NPTS2D_dum


lpj=0L

    readu, LUN[1], JMIN_IN ;(1:NLP_all)
if sw_debug eq 1 then print,' lpj=',lpj,' JMIN_IN(lpj)=',JMIN_IN;[0:5]
    readu, LUN[1], JMAX_IS ;(1:NLP_all)
if sw_debug eq 1 then print,'JMAX_IS(lpj)',JMAX_IS;[0:5],JMIN_IN[lpj+1]

if ( title_res eq 'td20120709' ) then  recalculate_jmin_max, JMIN_IN, JMAX_IS, NLP_ALL,sw_debug

;lpj=129L ;low
;lpj=34L ;dyn
;ipts=JMIN_IN(75)-1
ipts=JMAX_IS(75)-1

    dum = fltarr( NMP_all+1 ) ;!rad
    readu, LUN[1], dum ;( 1: NMP_all+1 ) !rad
if sw_debug eq 1 then print,'mlon[deg]',dum[mp]*180./!PI



    dum = Z_km
    readu, LUN[1], dum ;(     1:NPTS2D_dum) !meter
    Z_km     = dum * 1.0E-3
    if sw_debug eq 1 then print,lpj,ipts,' z_km',z_km[ipts]


if plot_type eq 5 then begin
   Re_km=6.3712E+03             ;km
   for i=0,203 do begin
      r=((z_km[i]+re_km)/re_km )
      print,'i=', i,' z_km=',z_km[i],' r=', r
;if r gt 2.5 then stop
   endfor ;i
;stop
endif ;projName


    dum  = mlat_deg
    readu, LUN[1], dum ;(    1:NPTS2D_dum) !rad
    mlat_deg = ( !PI*0.50 - dum ) * 180.0 / !PI
    if sw_debug eq 1 then print,'mlat_deg',mlat_deg[ipts]

    dum=fltarr(NPTS2D_dum, NMP_all)
    readu, LUN[1], dum
    glat_deg = ( !PI*0.50 - dum ) * 180.0 / !PI
    if sw_debug eq 1 then print,mp,'GCOLAT-deg',90.-dum[ipts, mp]*180./!PI


    dum=fltarr(NPTS2D_dum, NMP_all)
    readu, LUN[1], dum
    glon_deg = ( dum ) * 180.0 / !PI
    if sw_debug eq 1 then print,'GLON_rad-deg',dum[ipts, mp]*180./!PI  


;if sw_debug eq 1 then begin
   glatx= 42.8975
   glonx= 294.9765
   print, 'finding glatx=', glatx,'glonx=',glonx
   mp0=2
   lp0=31
   dmp=1
   dlp=1
   dlat=0.3
   dlon=4.7
   for mp=mp0-dmp, mp0+dmp,1 do begin
      if mp le 0 then mp1=mp+nmp0 else mp1=mp
      for lp=lp0-dlp, lp0+dlp,1 do begin

         if glatx ge 0. then $
            i=jmin_in[lp]-1L $ ;NH
         else $
            i=jmax_is[lp]-1L ;SH


print,mp1,lp, glat_deg[i,mp1],glon_deg[i,mp1]
         
         if ( (glat_deg[i,mp1]-dlat) lt glatx and glatx lt (glat_deg[i,mp1]+dlat) ) then begin 
            if ( (glon_deg[i,mp1]-dlon) lt glonx and glonx lt (glon_deg[i,mp1]+dlon) ) then begin 
               print,i,'mp1=',mp1,' lp=',lp,' glat=',glat_deg[i,mp1],' dlat=',(glat_deg[i,mp1]-glatx),' glon=',glon_deg[i,mp1],' dlon=',(glon_deg[i,mp1]-glonx),' z_km',z_km[i]
               print,' FLDIM=',(jmax_is[lp]-jmin_in[lp]+1)
            endif               ;glon
         endif                  ;glat
      endfor                    ;lp
   endfor                       ;mp
;STOP
   


   Re=6.3712E+03                ;km
   r_ref=Re+90.                 ;km for reference ht for rcm

   for lp=30,33 do begin
      i=jmin_in[lp]-1L

;assuming that latitude grid is identical for all LT sectors

      midpoint = JMIN_IN[lp]-1 + ( JMAX_IS(lp) - JMIN_IN(lp) )/2
      theta = (90.- mlat_deg[i])  *!PI / 180. ;[deg]-->[rad]
      sinthet = SIN( theta ) 
      lval    = r_ref / ( Re * sinthet * sinthet )

print,i,' lp=',lp,' mlat=',mlat_deg[i],' L-value=',lval;,(jmin_in(lp)-1),(jmax_is(lp)-1)
;print,midpoint,mlat_deg[midpoint]
   endfor
   ;stop

;endif ;if sw_debug eq 1 then begin


if ( sw_save_grid eq 1 ) then  begin
  print,'saving grid to a file=',filename_grid_sav
  save, filename=filename_grid_sav
  print,'saving grid finished'
endif



END ;PRO read_grid
