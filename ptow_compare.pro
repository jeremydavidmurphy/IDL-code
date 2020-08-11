PRO ptow_compare,file1, file2

readcol,file1,silent=1,f='d,d,d,d,d',c0,c1,c2,c3,c4

n0 = n_elements(c0)
wave1 = dblarr(2048,n0)

for k=0,n0-1 do begin
    for j=0,2047 do begin
        l = double(j)
        wave1[j,k] = c0[k] + c1[k]*l + c2[k]*l*l + c3[k]*l*l*l + c4[k]*l*l*l*l
    endfor
endfor

readcol,file2,silent=1,f='d,d,d,d,d',c0,c1,c2,c3,c4

n0 = n_elements(c0)
wave2 = dblarr(2048,n0)

for k=0,n0-1 do begin
    for j=0,2047 do begin
        l = double(j)
        wave2[j,k] = c0[k] + c1[k]*l + c2[k]*l*l + c3[k]*l*l*l + c4[k]*l*l*l*l
    endfor
endfor

speed = 299792.458

set_plot,'ps'
device,file='ptow_WAVEcomparison.ps',/portrait,/color
loadct,0,/silent
plot,wave1[*,0],wave1[*,0]-wave2[*,0],xrange=[3530,5900],/xstyle,$
  yrange=[-3,1.5],ytitle='Delta Wavelength (A)',$
  xtitle='Wavelength (A)',ystyle=1,xthick=2,ythick=2,charthick=2
;plot,wave1[*,0]/wave2[*,0],xrange=[0,2047],/xstyle,yrange=[0.9992,1.00002]
xyouts,0.20,0.3,file1,charsize=1.3,charthick=3,/normal
xyouts,0.20,0.25,'against',charsize=1.3,charthick=3,/normal
xyouts,0.20,0.2,file2,charsize=1.3,charthick=3,/normal
loadct,33
for j=0,n0-1 do oplot,wave1[*,j],wave1[*,j]-wave2[*,j],color=j
;for j=0,n0-1 do oplot,wave1[*,j]/wave2[*,j],color=j

device,/close_file

device,file='ptow_VELcomparison.ps',/portrait,/color
loadct,0,/silent
plot,wave1[*,0],((wave1[*,0]-wave2[*,0])/wave1[*,0])*speed,xrange=[3530,5900],/xstyle,$
  ytitle='Delta Velocity (km/sec)',yrange=[-40,40],$
  xtitle='Wavelength (A)',ystyle=1,xthick=2,ythick=2,charthick=2
;plot,wave1[*,0]/wave2[*,0],xrange=[0,2047],/xstyle,yrange=[0.9992,1.00002]
xyouts,0.20,0.3,file1,charsize=1.3,charthick=3,/normal
xyouts,0.20,0.25,'against',charsize=1.3,charthick=3,/normal
xyouts,0.20,0.2,file2,charsize=1.3,charthick=3,/normal
loadct,33
for j=0,n0-1 do oplot,wave1[*,0],((wave1[*,j]-wave2[*,j])/wave1[*,j])*speed,color=j
;for j=0,n0-1 do oplot,wave1[*,j]/wave2[*,j],color=j

device,/close_file
set_plot,'x'

stop
END
