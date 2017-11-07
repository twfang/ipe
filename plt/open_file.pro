;opening files
pro open_file,  input_DIR, LUN,version,input_flnm $
                ,sw_3DJ,sw_hr, sw_lun,title_res $
                ,sw_debug $
                , TimeStamp4Plot, n_read $
                , input_DIR0,luntmp7, luntmp3, luntmp8, luntmp9, luntmp10, luntmp11,sw_wam_ipe

  size_results=SIZE(input_DIR)
  n_file=size_results[1]
  if ( sw_debug eq 1 ) then  print,'n_file=',n_file

;flnmExt=getenv('flnmExt')

  if(sw_lun[0] eq 1 ) then  input_flnm[0]= $
     'ut_rec'+getenv('flnmExt') ;flnmExt
;'ut_recAll.ascii';concatinated
  if(sw_lun[1] eq 1 ) then  input_flnm[1]='plasma_grid.'+title_res
  if(sw_lun[2] eq 1 ) then  input_flnm[2]='plasma00' ;o+
  if(sw_lun[6] eq 1 ) then  input_flnm[6]='plasma01' ;h+
  if(sw_lun[8] eq 1 ) then  input_flnm[8]='plasma02' ;he+
  if(sw_lun[9] eq 1 ) then  input_flnm[9]='plasma03' ;N+
  if(sw_lun[10] eq 1 ) then  input_flnm[10]='plasma04' ;NO+
  if(sw_lun[11] eq 1 ) then  input_flnm[11]='plasma05' ;O2+
  if(sw_lun[12] eq 1 ) then  input_flnm[12]='plasma06' ;N2+
  if(sw_lun[13] eq 1 ) then  input_flnm[13]='plasma07' ;o+(2D)
  if(sw_lun[14] eq 1 ) then  input_flnm[14]='plasma08' ;o+(2P)
  
  if(sw_lun[3] eq 1 ) then  input_flnm[3]='plasma09' ;Te
  if(sw_lun[7] eq 1 ) then  input_flnm[7]='plasma10' ;Ti

  if(sw_lun[4] eq 1 ) then  input_flnm[4]='plasma12'  ;Vo+
  if(sw_lun[5] eq 1 ) then  input_flnm[5]='plasma16'  ;VExBup
  if(sw_lun[15] eq 1 ) then  input_flnm[15]='plasma17' ;VExBe
  if(sw_lun[16] eq 1 ) then  input_flnm[16]='plasma18' ;VExBth
  if(sw_lun[17] eq 1 ) then  input_flnm[17]='fort.2013' ;sunlons
;#nm20161128:plasma17 should be modified to plasma19 later...
  if(sw_lun[18] eq 1 ) then  input_flnm[18]='plasma17' ;sza
  if(sw_lun[20] eq 1 ) then  begin
     if n_read gt 0 then  free_lun, lun[20]
;     input_flnm[20]='ipe_grid_plasma_params.'+TimeStamp4Plot
;     ;george's new io  ;tmp20171025
     input_flnm[20]='ipe_grid_plasma_params';.'+TimeStamp4Plot  ;george's new io
     if ( sw_debug eq 1 ) then  print,'input_flnm[20]=', input_flnm[20]
  endif ;sw_lun20

;if ( sw_hr  eq 1 ) then begin
;for jth=1,7 do begin
;  input_flnm[jth+5]='fort.500'+STRTRIM( STRING( jth, FORMAT='(i1)'), 1) ;hrate(1) !!PGR 
;print,jth,(jth+5),input_flnm[jth+5]
;endfor
;endif
;if ( sw_3DJ eq 1 ) then input_flnm[6]='fort.4007' ;je3


  for i = 0, n_file-1  do begin

     if ( sw_lun[i] ne 1 ) then  CONTINUE

;dbg20121124      if ( i eq 0 ) or ( i gt 9 ) then $ 
     if ( i eq 0 OR i eq 17 )  then $ 
        openr, LUNi, input_DIR[i]+input_flnm[i], /GET_LUN $
     else $       
        openr, LUNi, input_DIR[i]+input_flnm[i], /GET_LUN $
               , /F77_UNFORMATTED

     if ( sw_debug eq 1 ) then  print,i,'luni=',luni
     LUN[i]=LUNi

     if ( i le 1 ) OR ( i eq 20 ) then       print,'opening file:', i,LUN[i],input_DIR[i]+input_flnm[i]
      
  endfor                        ;i = 0, n_file-1  do begin


  sw_read_wind = getenv('sw_read_wind')
  sw_version_io = getenv('sw_version_io')
  if ( sw_read_wind eq 1 ) then begin

     dirtmp7=input_DIR0
     print,'neutral DIR=',dirtmp7

     if sw_version_io eq 0 then $
        flnmtmp3 = dirtmp7+'ipe_grid_neut_wind' $ ;wind_out'       ;wam-ipe ;sw_version_io=0
     else if sw_version_io eq 1 then $
        flnmtmp3 = dirtmp7+'ipe_grid_neutral_params.'+TimeStamp4Plot ;sw_version_io=1 

;wind
     if n_read gt 0 then  free_lun, luntmp3
     openr, luntmp3, flnmtmp3, /GET_LUN $
            , /F77_UNFORMATTED
     print,'opening wind file=',flnmtmp3


     if sw_version_io eq 0 then begin
        luntmp7=201                              ;ut
        luntmp8=203                              ;tn
        luntmp9=204                              ;on
        luntmp10=205                             ;n2
        luntmp11=206                             ;o2
        flnmtmp7=dirtmp7+'ipe_grid_neut_params_ut' ;ut_out4wind' ;wam-ipe      
        flnmtmp8=dirtmp7+'ipe_grid_neut_temp'      ;tn_out'      ;tn
        flnmtmp9=dirtmp7+'ipe_grid_neut_O_den'     ;on_out'      ;on
        flnmtmp10=dirtmp7+'n2n_out'                ;n2
        flnmtmp11=dirtmp7+'o2n_out'                ;o2
        
      ;ut
        openr, luntmp7, flnmtmp7, /GET_LUN
        print,'opening ',flnmtmp7

;Tn
      openr, luntmp8, flnmtmp8, /GET_LUN $
            , /F77_UNFORMATTED
      print,'opening ',flnmtmp8
;on
      openr, luntmp9, flnmtmp9, /GET_LUN $
           , /F77_UNFORMATTED
      print,luntmp9,' opening ',flnmtmp9
;n2
;   openr, luntmp10, flnmtmp10, /GET_LUN $
;          , /F77_UNFORMATTED
;   print,'opening ',flnmtmp10
;o2
;   openr, luntmp11, flnmtmp11, /GET_LUN $
;          , /F77_UNFORMATTED
;   print,'opening ',flnmtmp11
     endif                      ;sw_version_io
  endif                         ;( sw_read_wind eq 1 ) then begin


end                             ;pro open_file
