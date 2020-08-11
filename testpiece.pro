pro testpiece

restore,'testpiece.idl'

medvalues = fltarr(n2,nfibers)
for j=0,n2-1 do begin
    for k=0,nfibers-1 do medvalues[j,k] = median(values[j,4,*,k])
endfor

window,0,retain=2
device,decomposed=0
loadct,0
plot,wavelength,values[*,4,0,0],yrange=[0.6,1.0],xrange=[330,620],/xstyle,/ystyle,psym=2,$
  /nodata,title='Fiber Transmission',xtitle='Wavelength (nm)',ytitle='% Transmission'
loadct,4

for j=0,nfibers-1 do begin
    for k=0,nstep-1 do oplot,wavelength,values[*,4,k,j],psym=sym(symbols[k]),color=colors[j]
    oplot,wavelength,medvalues[*,j],color=colors[j],psym=sym(12),symsize=2.5,thick=2
    fiber = strsplit(fiblist[j],'.',/extract)
    fiber = fiber[0]
    xyouts,0.75,0.5-(0.05*j),fiber,color=colors[j],/normal,charsize=1.5
endfor

stop
end
