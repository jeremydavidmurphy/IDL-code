FUNCTION resmapF, wave, IRarray, smth, med, name

; This code is used to generate a smooth map of the instrumental
; resolution. This is an altered version of resmap.pro. This is now
; a function that get's called once you have a wavelength array and IR
; array. It does the smoothing, then generates pretty 3D plots.

; It uses the FILTER_IMAGE routine.

; WAVE: The wavelength array
; IRARRAY: The IR array to be smoothed.
; SMTH: Size of the smoothing box.
; MED: Size of the median box.
;***********************************************************************
;smth = 3 ;the size of the smoothing box (MUST be odd). no higher than 5
;med = 3 ;the size of the median filter box (SHOULD be odd). keep below 5
;***********************************************************************

s = size(IRarray)
nn0 = s[1]
nn1 = s[2]
fwhm = IRarray

Pwave = fltarr(nn0,nn1)
for jj=0,nn0-1 do Pwave[jj,*] = replicate(wave[jj],nn1)

Py = findgen(nn1) + 1.0

izero = where(IRarray lt 0.0,countzero)
print,countzero
if (countzero gt 0) then IRarray[izero] = !Values.F_NAN
sfwhm = filter_image(IRarray,smooth=smth,median=med,/all_pixels)

zup = max(sfwhm)+0.1
zdn = min(sfwhm)-0.1
if (countzero gt 0) then sfwhm[izero] = -1.0

window,/free,retain=2,xsize=800,ysize=600
device,decomposed=0
loadct,0,/silent
shade_surf,IRarray,Pwave,Py,charsize=3.0,title='Original FWHM values',$
  zrange=[zdn,zup],/zstyle,ytitle='Fiber',ztitle='FWHM (A)',xtitle='Wavelength (A)'
;loadct,33,/silent
;shade_surf,IRarray,Pwave,Py,charsize=3.0,title='Original FWHM values',$
;  zrange=[zdn,zup],/zstyle,ytitle='Fiber',ztitle='FWHM (A)',xtitle='Wavelength (A)'
window,0,retain=2,xsize=800,ysize=600
device,decomposed=0
loadct,0,/silent
shade_surf,sfwhm,Pwave,Py,charsize=3.0,title='Smoothed FWHM values',$
  zrange=[zdn,zup],/zstyle,ytitle='Fiber',ztitle='FWHM (A)'
;loadct,33,/silent
;shade_surf,sfwhm,Pwave,Py,charsize=3.0,title='Smoothed FWHM values',$
;  zrange=[zdn,zup],/zstyle,ytitle='Fiber',ztitle='FWHM (A)'

form1 = '(a3,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x,f9.4,2x)'
form2 = '(i3,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6,2x)'
;form1 = '(a3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x)'
;form2 = '(i3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x,f9.3,2x)'
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
while !d.window ne -1 do wdelete, !d.window

out = 0.0
return,out
stop
END
