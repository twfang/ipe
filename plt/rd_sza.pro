;20161204underconstruction, test debug purpose
pro rd_sza

TEST='mpi20160330v2'
rundir='1480638394_ipe_theia_intel_parallel2_80'
rpath=$
'/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/runs/tmp20150730/trunk/run/'+rundir+'/'
;'~/stmp2/'+TEST+'/run/'+rundir+'/'

ut00=864000
utStart=ut00
utStop=utStart
NPTS2D=31287L
NMP=80L
restore, filename='~/ipeg/plt/plasma_grid.2xdyn.sav'

;o+
openr,lun00,rpath+'plasma00',/get_lun, /F77_UNFORMATTED
;sza
openr,lun01,rpath+'plasma17',/get_lun, /F77_UNFORMATTED
;ut_rec
openr,lun0, rpath+'ut_rec',/get_lun

dum=fltarr(NPTS2D,NMP) ;Te
dum1=fltarr(NPTS2D,NMP) ;sza
ut = 0L
record_number=0L
while ( eof(LUN00) eq 0 ) do begin
   readu, lun00,dum ;o+
   readu, lun01,dum1 ;sza
for  mp=0,nmp-1 do print,'sza=',mp,dum1[0,mp]*180./!PI
   readf, lun0,record_number, ut
   print,'rec#',record_number,' ut=', ut
;   ut = ut + dt ;[sec]
   if ut gt utStop then BREAK 
endwhile ;( eof(LUN00) eq 0 ) do begin

print,'rd_sza finished'

end ;pro rd_sza
