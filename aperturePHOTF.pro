FUNCTION aperturePHOTF,array,centX,centY,R1,R2,winscale

; This function conducts aperture photometry.

;ARRAY: The array
;centX: The x-position of the best guess for the center (in pixel
;coords, NOT IDL index)
;centY: The y-position of the best guess for the center
;R1: The inner radius, in pixels
;R2: The outer radius, in pixels

plottingF = 'on'

array = float(array)
centX = round(centX - 1.0)
centY = round(centY - 1.0)

chunk = array[centX-R2:centX+R2,centY-R2:centY+R2]
nn1 = n_elements(chunk[*,0])
nn2 = n_elements(chunk[0,*])

xx2 = indgen(R2)+1
xx1 = reverse(indgen(R2+1))
xline = [xx1,xx2]
yy2 = indgen(R2)+1
yy1 = reverse(indgen(R2+1))
yline = transpose([yy1,yy2])
Xarray = fltarr(nn1,nn2) & Yarray = Xarray

FOR i=0,nn1-1 DO Xarray(*,i) = xline
FOR i=0,nn2-1 DO Yarray(i,*) = yline
RADarr = sqrt(Xarray^2 + Yarray^2)

ispot = where(RADarr le R1)
iap = where(RADarr gt R1 and RADarr le R2)

if (plottingF eq 'on') then begin
    window,0,retain=2,xsize=nn1*winscale,ysize=nn2*winscale
    device,decomposed=0
    loadct,27,/silent
    tchunk = chunk
    tchunk[iap] = 0.0
    tvimage,chunk
    window,2,retain=2,xsize=nn1*winscale,ysize=nn2*winscale
    device,decomposed=0
    loadct,27,/silent
    tvimage,tchunk
endif

out1 = total(chunk[ispot])
out2 = mean(chunk[iap])
out3 = total(chunk[ispot]-out2)
out = [out1,out2,out3]
print,'The total pre-subtracted aperture flux is '+strn(out1)
print,'The mean annulus value is '+strn(out2)
print,'The total post-subtracted aperture flux is '+strn(out3)

return,out

END
