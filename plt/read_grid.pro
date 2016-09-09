;20121025: make sure to prepare the grid file (binary) in
;ipe/plt/grd/plasma_grid.XXX
pro read_grid,LUN,JMIN_IN,JMAX_IS,Z_km,mlat_deg,sw_debug,glat_deg,glon_deg,title_res


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

     mp=41-1L

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
;for i=0,nmp_all-1 do print,'mlon[deg]',(i+1),dum[i]*180./!PI
;for i=40,50 do print,'mlon[deg]',(i+1),dum[i]*180./!PI


    dum = Z_km
    readu, LUN[1], dum ;(     1:NPTS2D_dum) !meter
    Z_km     = dum * 1.0E-3
    if sw_debug eq 1 then print,lpj,ipts,' z_km',z_km[ipts]

;dbg20160811
;for i=0,50 do print, i,z_km[i]
;stop

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


glatx=42.6195
glonx=288.50827
print, 'glatx', glatx,'glonx',glonx
mp0=0
lp0=30;22;57;22
dmp=2;20
dlp=2;20
dlat=2.
dlon=3.
;for mp=mp0-dmp, mp0+dmp,1 do begin
;   for lp=lp0-dlp, lp0+dlp,1 do begin
;      i=jmin_in[lp]+35+15-10

;      if ( (glat_deg[i,mp]-dlat) lt glatx and glatx lt (glat_deg[i,mp]+dlat) ) then begin 
;         if ( (glon_deg[i,mp]-dlon) lt glonx and glonx lt (glon_deg[i,mp]+dlon) ) then begin 
;            print,i,'mp=',mp,' lp=',lp,' glat=',glat_deg[i,mp],(glat_deg[i,mp]-glatx),' glon=',glon_deg[i,mp],(glon_deg[i,mp]-glonx),' z_km',z_km[i]
;            print,jmin_in[lp],jmax_is[lp],' FLDIM=',(jmax_is[lp]-jmin_in[lp]+1)
;         endif                  ;glon
;      endif                     ;glat
;   endfor                       ;lp
;endfor                          ;mp
;STOP

;for lp=30,50 do begin
for lp=6-1,8-1 do begin
i=jmin_in[lp]
print,i,(lp+1),mlat_deg[i]
endfor
;stop

if ( sw_save_grid eq 1 ) then  begin
  print,'saving grid to a file=',filename_grid_sav
  save, filename=filename_grid_sav
  print,'saving grid finished'
endif



END ;PRO read_grid
