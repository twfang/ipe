

; Calculate TEC for altitude bands for a specified
; longitude and input times. 
; Altitude bands: 90 to 400 km, 400 to 800 km, 800 to 1000 km

pro tec_alt_bands

;Inputs
;1) root: Path to to directory containing ion density profile files
root='/ipe/smillholland/SSW_plasma_fixed/'   
;2) plot_root: Path to directory where output plot will be placed
plot_root='/ipe/smillholland/SSW_plasma_fixed/'
;3) times:Times to plot (LT hours) Local time is based on longitude
; indicated below.
times=[10,15]
;4) gLon: Geographic longitude to plot (West is negative)
gLon=-75
;5) time_diff: time difference such that LT=UT+time_diff based
; on longitude indicated above
time_diff=-5


; Open plasma files named stup00_fixed through.
;   stup08_fixed.
; 00: O+
; 01: H+
; 02: He+
; 03: N+
; 04: NO+
; 05: O2+
; 06: N2+
; 07: O+(2D)
; 08: O+(2P)

openr,lun0,root+'stup00_fixed',/GET_LUN
openr,lun1,root+'stup01_fixed',/GET_LUN
openr,lun2,root+'stup02_fixed',/GET_LUN
openr,lun3,root+'stup03_fixed',/GET_LUN
openr,lun4,root+'stup04_fixed',/GET_LUN
openr,lun5,root+'stup05_fixed',/GET_LUN
openr,lun6,root+'stup06_fixed',/GET_LUN
openr,lun7,root+'stup07_fixed',/GET_LUN
openr,lun8,root+'stup08_fixed',/GET_LUN

; Parameters of the ion density files
nT=97L ; number of times in the files
; the L says that the variables is a long integer
dT=900.0; increment between times (seconds)
startT=0.0 ;starting time (seconds)

nHt=183L ;number of height points
dHt = 5.0 ;increment between height points (km)
dHt_m = dHt*1000.0 ;increment between height points (m)
startHt = 90.0 ;starting height (km)
startHt_m = startHt*1000.0 ;starting height (m)

nLon=90L ;number of longitude points
dLon = 4.0 ; increment between longitude points (deg)
startLon = 0.0 ; starting longitude (deg)

nLat=91L ;number of latitude points
dLat = 2.0 ; increment bewteen latitude points (deg)
startLat = -90.0 ; starting latitude (deg)

; calculate height (km)
height=fltarr(nHt)
for i=0,nHt-1 do begin
  height[i]=(dHt*i)+ startHt
endfor
; calculate height (m)
height_m = height*1000.0

; calculate geographic latitude 
gLat=fltarr(nLat)
for k=0,nLat-1 do begin
  gLat[k]=(dLat*k)+startLat
endfor
; Take the low and mid latitudes (between -46 and 46 deg lat)
start=(-46-startLat)/dLat
stop=(46-startLat)/dLat
gLat=gLat[start:stop]

; Initialize TEC altitude band arrays
num_times=n_elements(times)
; nLat columns and num_times rows
TEC_90_400 = fltarr(nLat,num_times) 
TEC_400_800 = fltarr(nLat,num_times)
TEC_800_1000=fltarr(nLat,num_times)

; Calculate TEC altitudes band for specified times and 
;    specified longitude and for each latitude
; convert longitude to 0 to 360 deg coordinates
if (gLon lt 0) then gLon=gLon+360
; find index of longitude
lo=round((gLon-startLon)/dLon)
time_index=0 ; Index for times to calculate and plot
times=times-time_diff  ; Convert to UT
; Make sure no times are greater than 24
for t=0,num_times-1 do begin
   if (times[t] ge 24) then times[t]=times[t]-24
endfor

for t=0,nT-1 do begin    ; Time loop
   print,t 
   ; Define size of density profile arrays
   ion_dens0=fltarr(nHt,nLat,nLon) 
   ion_dens1=fltarr(nHt,nLat,nLon)
   ion_dens2=fltarr(nHt,nLat,nLon)
   ion_dens3=fltarr(nHt,nLat,nLon)
   ion_dens4=fltarr(nHt,nLat,nLon)
   ion_dens5=fltarr(nHt,nLat,nLon)
   ion_dens6=fltarr(nHt,nLat,nLon)
   ion_dens7=fltarr(nHt,nLat,nLon)
   ion_dens8=fltarr(nHt,nLat,nLon)
   elec_dens=fltarr(nHt,nLat,nLon)
   ; Read density profile arrays
   readf,lun0,ion_dens0
   readf,lun1,ion_dens1
   readf,lun2,ion_dens2
   readf,lun3,ion_dens3
   readf,lun4,ion_dens4
   readf,lun5,ion_dens5
   readf,lun6,ion_dens6
   readf,lun7,ion_dens7
   readf,lun8,ion_dens8
   ; The arrays are summed to calculate electron 
   ;   density.
   elec_dens=ion_dens0+ion_dens1+ion_dens2+ion_dens3+ion_dens4 $
      +ion_dens5+ion_dens6+ion_dens7+ion_dens8
   time=(t*dT)/3600.0 ; Time in hours
   result=where(times eq time)
   if (result ne -1) then begin
   for la = 0, nLat-1 do begin       ; Latitude loop
      density = elec_dens[*,la,lo] 
      ; Calculate TEC altitude bands in units of TECU
      ; Find indices of heights
      result=min(abs(height-90),i_90)
      result=min(abs(height-400),i_400)
      result=min(abs(height-800),i_800)
      result=min(abs(height-1000),i_1000)
      ; Partition density
      dens_90_400=density[i_90:i_400]
      dens_400_800=density[i_400:i_800]
      dens_800_1000=density[i_800:i_1000]
      ; Partition height
      height_m_90_400=height_m[i_90:i_400]
      height_m_400_800=height_m[i_400:i_800]
      height_m_800_1000=height_m[i_800:i_1000]
      ; Calculate TEC 
      TEC_90_400_point = int_tabulated(height_m_90_400,$
         dens_90_400)/1E16
      TEC_400_800_point = int_tabulated(height_m_400_800,$
         dens_400_800)/1E16
      TEC_800_1000_point = int_tabulated(height_m_800_1000,$
         dens_800_1000)/1E16 
      TEC_90_400[la,time_index] = TEC_90_400_point
      TEC_400_800[la,time_index] = TEC_400_800_point
      TEC_800_1000[la,time_index] = TEC_800_1000_point
   endfor       ; Latitude loop
   
   ; Advance time index for specified times to calculate and plot 
   time_index=time_index+1
   print,'time index',time_index
   endif ; If statement which determined whether to calculate and plot
   if (time_index eq num_times) then break
endfor    ; Time loop


; Convert the bands to additive arrays
TEC_800_1000=TEC_90_400+TEC_400_800+TEC_800_1000
TEC_400_800=TEC_90_400+TEC_400_800


; Take the low and mid latitudes (between -46 and 46 deg lat)
start=(-46-startLat)/dLat
stop=(46-startLat)/dLat
if (num_times eq 1) then begin
   TEC_90_400=TEC_90_400[start:stop]
   TEC_400_800=TEC_400_800[start:stop]
   TEC_800_1000=TEC_800_1000[start:stop]
endif else begin
   TEC_90_400=TEC_90_400[start:stop,*]
   TEC_400_800=TEC_400_800[start:stop,*]
   TEC_800_1000=TEC_800_1000[start:stop,*]
endelse

; Close files
close,lun0
close,lun1
close,lun2
close,lun3
;close,lun4
;close,lun5
;close,lun6
;close,lun7
;close,lun8


; Create TEC altitude band plots for each time

; Plotting commands 
;Enable colors
device,decomposed=0,retain=2 
window,/free,/pixmap,colors=256&wdelete,!d.window
; Create contour colorbars for NmF2, HmF2, and TEC
set_plot, 'PS'
aspect_ratio=1
xsize=15
ysize=xsize/aspect_ratio
!p.font=0
device,xsize=xsize,ysize=ysize
loadct,12

; Plot and save files

; Convert gLon back to -180 to 180 degree coordinates 
; for plotting
if (gLon gt 180) then gLon=gLon-360

; Convert times back to LT for plotting
times=times+time_diff
; Make sure no times are less than 0
for t=0,num_times-1 do begin
   if (times[t] le 0) then times[t]=times[t]+24
endfor


; Begin plotting
for t=0,num_times-1 do begin ; Plotting loop
   time=times[t]
   TEC_90_400_1=TEC_90_400[*,t]
   TEC_400_800_1=TEC_400_800[*,t]
   TEC_800_1000_1=TEC_800_1000[*,t]

   ; Calculate percentages
   ; Integrated curves
   int_90_400=int_tabulated(gLat,TEC_90_400_1)
   int_400_800=int_tabulated(gLat,TEC_400_800_1)
   int_800_1000=int_tabulated(gLat,TEC_800_1000_1)
   ; Integrated percentages
   per_int_90_400=(int_90_400/int_800_1000)*100
   per_int_400_800=((int_400_800-int_90_400)/int_800_1000)*100
   per_int_90_400=round(per_int_90_400*10)/10.0
   per_int_400_800=round(per_int_400_800*10)/10.0
   per_int_800_1000=100.0-per_int_90_400-per_int_400_800

   ; Non-integrated percentage of curve below 400 km
   per_90_400=(TEC_90_400_1/TEC_800_1000_1)*100


   ; Create plots
   device,filename=plot_root+'TEC_alt_bands_long='+strtrim(gLon,2)+$
   'time='+strtrim(time,2)+'LT.ps',$
   encapsulated=1,/helvetica,/color
   ; Plot altitude bands
   plot,gLat,TEC_800_1000_1,xtitle = 'Latitude (deg)',ytitle = $
      'TEC (TECU)',title = 'TEC Altitude Bands: LT='+strtrim(time,2)+$
      '(hours) Longitude='+strtrim(gLon,2)+' (deg)',$
      xthick = 2,ythick = 2,xrange=[gLat[0],gLat[-1]],$
      yrange=[0,35],position = [0.2,0.4,0.76,0.9],$
      xstyle=1,ystyle=8,/nodata
   pxval=[gLat[0],gLat,gLat[-1]]
   ; Fill in with red, green, and blue
   polyfill,pxval,[0,TEC_800_1000_1,0],color=104
   polyfill,pxval,[0,TEC_400_800_1,0],color=40
   polyfill,pxval,[0,TEC_90_400_1,0],color=200
   

   ; Create new axis for plot of percentage of curve below 400 km
   axis,yaxis=1,yrange=[0,100],ytitle=$
      'Percentage of TEC Below 400 km',$
      xthick=2,ythick=2,ticklen=-0.02,color=120,/save
   ; Plot percentage 
   oplot,gLat,per_90_400,color=120
   

   ; Create legend
   levels_contour=indgen(17)*16
   arr=fltarr(2,90)
   arr[*,0:29]=arr[*,0:29]+200+16
   arr[*,30:59]=arr[*,30:59]+40+16
   arr[*,60:89]=arr[*,60:89]+104+16
   contour,arr,$
      levels=levels_contour,position = [0.89,0.5,0.92,0.8],/noerase,$
      /fill,xrange=[0.0,2.0],yrange=[0,90],xstyle=4,ystyle=4
   xyouts,1.5,15,'90 to 400 km'
   xyouts,1.5,45,'400 to 800 km'
   xyouts,1.5,75,'800 to 1000 km'
   ; Print percentages
   xyouts,1.5,5,$
      strmid(strtrim(per_int_90_400,2),0,4)+'%'
   xyouts,1.5,35,$
      strmid(strtrim(per_int_400_800,2),0,4)+'%'
   xyouts,1.5,65,$
      strmid(strtrim(per_int_800_1000,2),0,4)+'%'

endfor ; Plotting loop
device,/close 

print, 'Procedure finished successfully!'
end ;pro read_ion_densities
