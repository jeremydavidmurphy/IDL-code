PRO input_check, list

; This routine is used to confirm the input f-ratio of the test
; bench. It accepts a list of fits files of the following form:
; camera position     diameter

cstep = 2.0 ; the camera step size, in mm
psize = 0.0120000; the pixel size, in mm
thresh = 2000.0
winsize=4.0
;*********************************************************************

readcol,list,f='f,f',dis,dia
dia = dia - dia[0]
dia = dia * psize * 2.0
n0 = n_elements(dis)
window,0,retain=2

line = linfit(dis,dia)
xarray = (1.0+findgen(100))*(max(dis)/100.0)
xarray = [0,xarray]
yarray = fltarr(101)

for j=0,n_elements(xarray)-1 do yarray[j] = line[0] + line[1]*xarray[j]

fout = xarray / yarray

loadct,0
plot,dis,dia,psym=1,yrange=[0,6],/ys
loadct,4
oplot,xarray,yarray,color=150
xyouts,0.09,0.8,'f/in = ',charsize=4,charthick=2,/normal
xyouts,0.12,0.8,mean(fout[50:*]),charsize=4,charthick=2,/normal

set_plot,'ps'
device,file=list+'.ps',/color
loadct,0
plot,dis,dia,psym=1,yrange=[0,6],/ys,xthick=3,ythick=2,charthick=2,$
  xtitle='Travel (mm)',ytitle='Cone Size (mm)',charsize=1.2
loadct,4
oplot,xarray,yarray,color=150,thick=2
xyouts,0.19,0.8,'f/in = ',charsize=2,charthick=2,/normal
xyouts,0.22,0.8,mean(fout[50:*]),charsize=2,charthick=2,/normal
device,/close
set_plot,'x'

stop
END
