PRO resmap, wavearray, IRarray

; This code is used to generate a smooth map of the instrumental
; resolution.  based on fitted gaussians to an extracted arc frame. The
; returned values are the text files, month_night_FWHMraw.txt and
; month_night_FWHMsmoothed.txt values. This is a per fiber, per
; wavelength estimation (both raw and smoothed) of the FWHM of the
; gaussian fit to the arcs for that night.

; It uses the FILTER_IMAGE routine.

; NAME: This is the name AND NIGHT of the run. When fed into the primary
; workhorse of this code (rescheckF) the extracted arc, ptow.dat and
; mask.dat file are located as below.
; (the '_' is subscripted onto the name.)
;        arc = 'arc'+name+'e.fits'
;        ptow = 'ptow'+name
;        mask = 'mask'+name

; COMPILE: rescheckF

;***********************************************************************
smth = 3 ;the size of the smoothing box (MUST be odd). no higher than 5
med = 3 ;the size of the median filter box (SHOULD be odd). keep below 5
;***********************************************************************

values = rescheckF('_'+name)

;wave = [4046.5539, 4077.8298, 4358.3253, 4678.1474, 4799.9080, 4916.0661,$
;         5085.8173, 5154.6616, 5460.7366, 5769.6000, 5790.6239]
wave = [4046.5539, 4077.8298, 4358.3253, 4678.1474, 4799.9080,$
         5085.8173, 5460.7366, 5769.6000, 5790.6239]

fwhm = values[*,*,0]
disp = values[*,*,1]

sfwhm = filter_image(fwhm,smooth=smth,median=med,/all_pixels)
sdisp = filter_image(disp,smooth=smth,median=med,/all_pixels)

zup = max(sfwhm)+0.1
zdn = min(sfwhm)-0.1
window,2,retain=2,xsize=800,ysize=600
device,decomposed=0
shade_surf,fwhm,charsize=3.0,title='Original FWHM values',$
  zrange=[zdn,zup],/zstyle,ytitle='Fiber',ztitle='FWHM (A)'
window,0,retain=2,xsize=800,ysize=600
device,decomposed=0
shade_surf,sfwhm,charsize=3.0,title='Smoothed FWHM values',$
  zrange=[zdn,zup],/zstyle,ytitle='Fiber',ztitle='FWHM (A)'

;form1 = '(f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x)'
;form2 = '(f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x)'
form1 = '(a3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x)'
form2 = '(i3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x)'
free_lun,5
openw,5,name+'_FWHMsmoothed.txt'
printf,5,'F#',wave,format=form1
for j=0,n_elements(sfwhm[0,*])-1 do printf,5,j+1,sfwhm[*,j],format=form2
free_lun,5
openw,5,name+'_FWHMraw.txt'
printf,5,'F#',wave,format=form1
for j=0,n_elements(fwhm[0,*])-1 do printf,5,j+1,fwhm[*,j],format=form2
free_lun,5

print,'Next ENTER deletes the plots...'
pause
wdelete
wdelete

stop
END
