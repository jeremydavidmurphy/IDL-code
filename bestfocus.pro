PRO bestfocus, list

; This routine reads in the fibernam_focus.txt files that come out of
; fiber_trans2 routine and returns which frame is the best focus.

; The fiber focus text frames contain 2 columns. The first is the median
; radius of the spot size. This is the median of 21 estimates up and down the
; height of the fiber cross-sectional profile. Although this should correlate
; with best focus, there are issues with taking this number directly
; (i.e. just finding the smallest fiber spot diameter and calling that best
; focus). This is due to the fact that the caluclation of the spot diameter is
; based on a percentage of the total amount of light in the fiber. Therefore,
; there's dependencies on total flux, evenness of illumination, weak
; fibers, etc. The second column is the standard deviation (SD) of the 21
; estimates. This allows for an exploration of the sharpness of the profile
; rather than the diameter. The idea is that a well-focused image will
; have a very-nearly top-hat fiber profile. The SD of the 21
; measurements of the spot size at different levels of encircled
; energy (EE) gives a measure of how steep or shallow that slope
; is. If the image is in perfect focus then the diameter at all EE
; will be the same and the SD will be zero. 

readcol,list,silent=1,f='a',files
n0 = n_elements(files)

for j=0,n0-1 do begin
   readcol,files[j],silent=1,f='f,f',a,b
   if (j eq 0) then array = fltarr(2,n_elements(a),n0)
   array[0,*,j] = a
   array[1,*,j] = b
endfor

window,2,retain=2
device,decomposed=0
loadct,0,/silent
;plot,array[0,*,0],psym=-1,/nodata,yrange=[56,62],/ystyle,$
;     ytitle='Radius measurements (pixels)'
plot,array[1,*,0],psym=-1,/nodata,yrange=[3,8],/ystyle,charsize=1.3,$
     ytitle='Standard Deviation of Radius Measurements (pixels)'
loadct,33,/silent
for j=0,n0-1 do begin
;   oplot,array[0,*,j],psym=-1,color=j*20
   oplot,array[1,*,j],psym=-1,color=50+(j*20)
   print,files[j]
   pause
endfor

;best = median(array,dim=2)
;best = best[0,*]
best = median(array,dim=2)
best = best[1,*]

ibest = where(best eq min(best))
print,'The best focus is found at frame '+files[ibest]

stop
END






