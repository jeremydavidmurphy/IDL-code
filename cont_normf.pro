; This routine takes a spectra and runs a boxcar over it. It then
; returns both the smoothed version, and the continuum divided
; version.

FUNCTION cont_normf, data, BOX=box

;***************************************************************
plottingF = 'off'
;plottingF = 'on'
;***************************************************************

if (n_elements(box) eq 0) then box = 150

dsmtd = smooth(data,box,/edge_truncate,/nan)
sdata = data / dsmtd

ihigh = where(sdata gt (median(sdata) * 2.0),count)
if (count ne 0) then begin
    data[ihigh] = !Values.F_NAN
    dsmtd = smooth(data,box,/edge_truncate,/nan)
    sdata = data / dsmtd
endif

if (plottingF eq 'on') then begin
    window,4,retain=2
    device,decomposed=0
    loadct,0,/silent
    plot,data,xrange=[30,2000],/xstyle,/ynozero
    loadct,4,/silent
    oplot,dsmtd,color=150
    pause
    loadct,0,/silent
    wdelete,4
endif

out = [[sdata],[dsmtd]]
return,out

END
