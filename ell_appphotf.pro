FUNCTION ell_appphotf,array,centX,centY,pa,major,minor,R1f,R2f,TYPE=type,WIN=winscale,PLOT=plot

; This routine is used by fiber_trans2.pro to conduct elliptical aperture
; photometry. It estimates the background in an annulus, using either the
; median or mean of the counts in the aperture. It then subtracts this value
; from every pixel within the inner radius, then totals the counts within the
; radius.

;ARRAY: The array
;centX: The x-position of the best guess for the center (in pixel
;coords, NOT IDL index)
;centY: The y-position of the best guess for the center
;pa: The position angle of the ellipse. Measured CCW from the Y-axis
;major: The ellipse major axis
;minor: The ellipse minor axis
;R1f: The amount of radius to add to the major axis to enlarge the included region.
;R2f: The amount to add to the major radius to estimate sky (the ellipse gt R1f
;and lt R2f forms the sky estimate
;winscale: The scaling for the window size
;type: set to either 'median' or 'mean'. If the keyword isn't applied,
;the mean is taken
;plot: set this to ANYTHING and the plotting will happen

if (n_elements(plot) eq 0) then plottingF = 'on' else plottingF = plot

if (n_elements(type) eq 0) then type = 'mean' ;will use the mean, if not instructed to use the median
if (n_elements(winscale) eq 0) then winscale = 10.0

ratio = float(major)/float(minor)
array = float(array)
centX = round(centX - 1.0)
centY = round(centY - 1.0)
majorC = round(major)

chunk = array[centX-majorC-R2f:centX+majorC+R2f,centY-majorC-R2f:centY+majorC+R2f]
nx = n_elements(chunk[*,0])
ny = n_elements(chunk[0,*])

;**************************************************************************
; The ellipse array is created. This is akin to the radius array for the
; circular aperture photometry routine. The conversions and code follows the
; dist_ellipse.pro routine.

ang = pa / !RADEG ;convert position angle from degrees to radians
cosang = cos(ang)
sinang = sin(ang)

xf = findgen(nx) - (nx/2.0)
yf = findgen(ny) - (ny/2.0)
ellarr = fltarr(nx,ny,/nozero)
xcosang = xf*cosang ;rotate pixels to match ellipse orientation
xsinang = xf*sinang

for ii = 0,ny-1 do begin
   xtemp =  xcosang + yf[ii]*sinang
   ytemp = -xsinang + yf[ii]*cosang
   ellarr[0,ii] = sqrt( (xtemp*ratio)^2 + ytemp^2 )
endfor
;**************************************************************************

i1 = where(ellarr le major+R1f)
i2 = where(ellarr gt major+R1f and ellarr le major+r2)
ellarr2 = ellarr
ellarr2[i2] = 0.0
pchunk = chunk
pchunk[i2] = 0.0

if (plottingF eq 'on') then begin
    window,3,retain=2,xsize=nx*winscale,ysize=ny*winscale
    device,decomposed=0
    loadct,37,/silent
    tvimage,bytscl(pchunk,top=!d.table_size-3)
endif

;Now, the background is subtracted...
if (type eq 'median') then M1 = median(chunk[i2],/even)
if (type eq 'mean') then M1 = mean(chunk[i2])
chunkout = chunk - M1
out = [total(chunkout),M1]

return,out ;the total, and the background per pixel are subtracted.

END
