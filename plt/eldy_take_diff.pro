pro eldy_take_diff

RDIR00='/scratch3/NCEPDEV/swpc/scrub/Naomi.Maruyama/fig/efield/20171025/raw_high_lat/'
RDIR0=RDIR00+'JQIPEr420/'
RDIR1=RDIR00+'IPEOptimization/'
jmax=1L
for j=0,jmax-1 do begin
   restore,RDIR0+'eldynUTsec0.serial.sav'
help
   print,'(0) utime=', utime

   if j eq 0 then $
      var0=poten $
   else if j eq 1 then $
      var0=ed1130 $
   else if j eq 2 then $
      var0=ed2130 $
   else if j eq 3 then $
     var0=ed190 $
   else if j eq 4 then $
      var0=ed290
   
   max0=MAX(var0)
   min0=Min(var0)
   ave0=(max0+min0)*0.5
   print,j,'(0)poten:MAX=',max0,min0
   

   restore,RDIR1+'eldynUTsec0.parallel_1.sav'
;help
   print,'(1) utime=', utime
   

   if j eq 0 then $
      var1=poten $
   else if j eq 1 then $
      var1=ed1130 $
   else if j eq 2 then $
      var1=ed2130 $
   else if j eq 3 then $
      var1=ed190 $
   else if j eq 4 then $
      var1=ed290

   max1=MAX(var1)
   min1=Min(var1)
   print,j,'(1)poten:MAX=',MAX1,MIN1

   diff=var1-var0
   print,j,'diff=',max(diff),min(diff)
   print,j,'rat=',max(diff)/max0, min(diff)/max0
   
endfor                          ;j

print, "eldy_take_diff finished!"
end;pro eldy_take_diff
