;date 20160827
;purpose: how to asign gravity waves perturbations to wind in IPE???
pro plt_sine_wave

imax=7200L
y=fltarr(imax)
x=fltarr(imax)
A=20.
f=1./3600.
phi=0.
t=findgen(imax) ;[sec]
print,'t[sec]',t
for i=0,imax-1 do begin
   x[i]=2. * !PI * f * t[i] + phi
   y[i]= A * SIN(x[i])
endfor
;plot, x,y
plot, t/3600.,y

loadct, 39
phi=!PI*2.
for i=0,imax-1 do begin
   x[i]=2. * !PI * f * t[i] + phi
   y[i]= A * SIN(x[i])
endfor
;oplot, x,y, linestyle=2, color=240
oplot, t/3600.,y, linestyle=2, color=240
end;pro plt_sine_wave
