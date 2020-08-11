PRO astrom_compare, list
  
;this routine is used to compare the astrometry of the fiber positions
;for several repeat measurements. The LIST is a list of each of the
;text files that's produced by fiber_trans2. The first file in
;this list is assumed to be the fiducial, from which all other
;positions are referenced.

;I have calculated the plate scale for the reimager, from both VP2.7
;astrometry data and an VIRUS IFU head. The relative scales calculated
;are:
; From the VIRUS IFU: 4.354 microns/pixel
; From the VP2.7 IFU: 4.345 microns/pixel
; This difference could be simply because the VIRUS IFU was rastered
;around the FOV and so would have suffered from slightly larger
;distortion.

readcol,list,f='a',silent=1,files
n0 = n_elements(files)

for j=0,n0-1 do begin
    readcol,files[j],silent=1,f='i,d,d,d,x,x,x,x',skipline=1,a,b,c,d
    if (j eq 0) then begin
        n1 = n_elements(a)
        fibers = a
        Xarr = dblarr(n1,n0)
        Yarr = dblarr(n1,n0)
        Rarr = dblarr(n1,n0)
        Rfirst = d * 4.35 ;from pixels to microns
    endif
    Xarr[*,j] = b
    Yarr[*,j] = c
    Rarr[*,j] = Rfirst - (d * 4.35) ;from pixels to microns
endfor

cmulti = floor(255.0/n0)
window,0,retain=2
device,decomposed=0
loadct,0,/silent
plot,a,Rarr[*,0],psym=1,yrange=[-1.5,1.5],/ys
loadct,33,/silent
wait,0.3
for j=1,n0-1 do begin
    oplot,a,Rarr[*,j],psym=1,color=j*cmulti
    wait,0.2
endfor

set_plot,'ps'
device,file='astrom_compare.ps',/color
loadct,0,/silent
plot,a,Rarr[*,0],psym=1,yrange=[-1.5,1.5],/ys,xtitle='Fiber #',$
  ytitle='Delta radial distance, fiber-to-fiber center (microns)',$
  xthick=2,ythick=2,charthick=2,xrange=[0,225],/xs,/nodata
loadct,33,/silent

for j=1,n0-1 do begin
    oplot,a,Rarr[*,j],psym=1,color=j*cmulti,thick=2
endfor
device,/close
set_plot,'x'

stop
END
