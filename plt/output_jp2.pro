;20150811: read image from the current display and write it to the jpeg2000(jp2) file

PRO output_jp2, Filename_jp2


; readin the current WINDOW as image
thisImage= TVRD( $
;[X 0 [, Y 0 [, N x [, N y [, Channel]]]]] [, CHANNEL=value] $
; [,/ORDER] [ $
  TRUE=1  $ ;{1 | 2 | 3}]  ;110904!CAUTION! 2or3 did not work!!! 
;[, /WORDS] $
)

;check if the image was properly read
;20150811 image does not appear!!! i dont understand why???
;WINDOW,3,XSIZE=900,YSIZE=700
;tv, thisImage


; output to jp2 file
WRITE_JPEG2000, $
Filename_jp2 $
, thisImage  $
;[, Red, Green, Blue] $
;note: default is 1
, N_LAYERS=1 ;$
;[, N_LEVELS=value] [, /ORDER] [, /REVERSIBLE]

; WRITE_PNG  $
; , Filename_png $
; , thisImage   ;[, R, G, B] 
;;[, /ORDER] [, /VERBOSE] [, TRANSPARENT=array]


 end ;PRO output_jp2
