pro read_fort167,LUN167,UT_hr,LT_hr,Z,SL,GL,BM,GR,SZA,O0,H0,N20,O20,HE0,N4S0,sw_debug, title_hemi


  ;get FLDIM
   size_result = size(Z)
if ( sw_debug eq 1 ) then   print, 'size=',size_result
   FLDIM = size_result[1]


   string_tmp1='U167, North, UT='
;NH
   readf, LUN167, string_tmp1, UT_hr, LT_hr, FORMAT='(A16,2F10.2)'
if ( sw_debug eq 1 ) then      print, string_tmp1, UT_hr, LT_hr


   string_tmp2='   JMIN   JMAX   CTIPDIM  INNO  IHEPLS  INPLS'
JMIN=0L
JMAX=0L
CTIPDIM=0L
INNO=0L
IHEPLS=0L
INPLS=0L
;NH
   readf, LUN167, string_tmp2
if ( sw_debug eq 1 ) then      print,   'NH',string_tmp2
   if ( title_hemi eq 'NH' ) then begin
     readf, LUN167, JMIN,  JMAX,  CTIPDIM, INNO, IHEPLS, INPLS, FORMAT='(22I7)'
     if ( sw_debug eq 1 ) then      print,   'NH',JMIN,  JMAX,  CTIPDIM, INNO, IHEPLS, INPLS
   endif

   string_tmp3='        PCO       DT          DTMIN       F107       F107A          FPAS        HPEQ      HEPRAT      COLFACX'

;blank line
;NH
   readf, LUN167, string_tmp3
if ( sw_debug eq 1 ) then      print, 'NH blank1!',string_tmp3


;NH
   readf, LUN167, string_tmp3
if ( sw_debug eq 1 ) then      print, 'NH',string_tmp3

   if ( title_hemi eq 'NH' ) then begin
     readf, LUN167,PCO,      DT,         DTMIN,      F107,      F107A,         FPAS,       HPEQ,     HEPRAT,     COLFACX, FORMAT='(22F12.3)'
     if ( sw_debug eq 1 ) then      print, 'NH',PCO,      DT,         DTMIN,      F107,      F107A,         FPAS,       HPEQ,     HEPRAT,     COLFACX
   endif


   string_tmp4='     Z      SL       GL      BM        GR       SZA        O         H       N2       O2       HE       N4S'

;blank line
;NH
   readf, LUN167, string_tmp4
if ( sw_debug eq 1 ) then      print, 'NH blank2!',string_tmp4



;NH
   readf, LUN167, string_tmp4
if ( sw_debug eq 1 ) then      print, 'NH',string_tmp4





   for j=1-1,(FLDIM/2)+1-1 DO begin
;NH
   readf, LUN167, Zj,SLj,GLj,BMj,GRj,SZAj,Oj,Hj,N2j,O2j,HEj,N4Sj , FORMAT='(F10.2,E14.7,21E9.2)'

  Z(J) =zj
 SL(J) =SLj
 GL(J) =GLj
 BM(J) =BMj
 GR(J) =GRj
SZA(J) =SZAj
  O0(J) =Oj
  H0(J) =Hj
 N20(J) =N2j
 O20(J) =O2j
 HE0(J) =HEj
N4S0(J) =N4Sj


;if ( sw_debug eq 1 ) then   if ( j eq 10 ) then  $
;print, j, Zj,SLj,GLj,BMj,GRj,SZAj,Oj,Hj,N2j,O2j,HEj,N4Sj , FORMAT='(i4,F10.2,22E9.2)'

   ENDFOR

END ;pro read_fort167
