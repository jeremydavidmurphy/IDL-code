PRO Prmsall, rmslist

; This routine is used to plot the rmsall.out.?? files generated when
; you run rfitlov and rfitlovc.

readcol,rmslist,silent=1,format='a',rmsfiles
n0 = n_elements(rmsfiles)

for l=0,n0-1 do begin
    readcol,rmsfiles[l],silent=1,format='f,x,x,x',rms
    n1 = n_elements(rms)
    if (l eq 0) then begin
        window,0,retain=2
        device,decomposed=0
        loadct,0
        plot,rms,psym=1,xrange=[-2,n1+2],xstyle=1,/nodata,title='RMS VALUES',$
          xtitle='Files',ytitle='RMS',charsize=1.2
        loadct,4
    endif
    oplot,rms,psym=-1,color=60+(l*50)
    xyouts,0.15,0.90-(l*0.04),rmsfiles[l],color=60+(l*50),/normal,charsize=1.2
endfor

print,'Take a look...'
print,'Next ENTER deletes the window...'
pause
wdelete

stop
END
