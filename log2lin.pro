FUNCTION log2lin, spectrain, logwave1, logwavedisp, linwave1, linwavedisp

; reads in a fits frame that's in log wavelength space and
; outputs a spectra that's linear wavelength space. The log
; wavelength is assumed to be of the form:
;
; lambda_x = logwave1 * 10^(logwavedisp * x)
;
; where logwave1 is the value of pixel position 1 (e.g. the 'x' array
; starts with zero)

n0 = n_elements(spectrain)

x = dindgen(n0)

wavelog = logwave1 * 10.0^(logwavedisp * x)
wavelin = dindgen(n0) * linwavedisp + linwave1

spectraout1 = interpol(spectrain,wavelog,wavelin)
spectraout2 = interpol(spectrain,wavelog,wavelin,/quadratic)
spectraout3 = interpol(spectrain,wavelog,wavelin,/lsquadratic)
spectraout4 = interpol(spectrain,wavelog,wavelin,/spline)

;window,0,retain=2
;device,decomposed=0
;loadct,0
;plot,wavelin,spectraout1,yrange=[0,2],xrange=[3600,5400]
;loadct,4
;pause
;oplot,wavelin,spectraout2,color=60
;pause
;oplot,wavelin,spectraout3,color=110
;pause
;oplot,wavelin,spectraout4,color=150
;pause
;loadct,0
;oplot,wavelog,spectrain

return,spectraout4

stop
END



