PRO input_check_radius, list

; This routine is used to confirm the input f-ratio of the test
; bench.

; LIST: The list of fits files.

cstep = 2.0 ; the camera step size, in mm
psize = 0.0120000; the pixel size, in mm
thresh = 4000.0

a1 = 3.14159*(606.515^2)
a2 = 3.14159*(262.822^2)
;a1 = 3.14159*(776.0^2)
;a2 = 3.14159*(330.6^2)
trim = a2/a1
;*********************************************************************

window,2,retain=2,xsize=600,ysize=600

readcol,list,f='a',files
n0 = n_elements(files)
radius = fltarr(n0)
steps = findgen(n0)*cstep
for j=0,n0-1 do begin
    file = files[j]
    frame = readfits(file,/silent)
    ihigh = where(frame gt thresh,nhigh)
    frame[ihigh] = 0.0
    tvscale,frame
    area1 = psize^2 * nhigh
    area1 = (1 + trim) * area1
    radius[j] = (sqrt(area1/3.14159))*2.0
    print,radius[j]
endfor

radius = radius - radius[0]
line = linfit(steps,radius)
xarray = (1.0+findgen(100))*(max(steps)/100.0)
xarray = [0,xarray]
yarray = fltarr(101)

for j=0,n_elements(xarray)-1 do yarray[j] = line[0] + line[1]*xarray[j]
fout = xarray / yarray

window,0,retain=2
device,decomposed=0
plot,steps,radius,psym=4
loadct,4
oplot,xarray,yarray,color=150
xyouts,0.09,0.8,'f/in = ',charsize=4,charthick=2,/normal
xyouts,0.12,0.8,mean(fout[50:*]),charsize=4,charthick=2,/normal

stop
END
