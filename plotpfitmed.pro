pro plotpfitmed, list

yupvel = 100
ydownvel = -200
yupsig = 600
ydownsig = 250

color = [60,75,90,105,120,135,150,165,180,195,210,225,240,255]
symbol = [1,2,4,5,6,7,1,2,4,5,6,7]
readcol,list,format='a',datalist

n1 = n_elements(datalist)
readcol,datalist[0],format='i,i',vel,sig
n2 = n_elements(vel)

velarr = intarr(n1,n2)
sigarr = intarr(n1,n2)

for j=0,n1-1 do begin
    readcol,datalist[j],format='i,i',vel,sig
    velarr[j,*] = vel
    sigarr[j,*] = sig
endfor

window,0,retain=2,xsize=900
device,decomposed=0

loadct,0
plot,velarr[0,*],/ynozero,/nodata,title='VELOCITY',$
  yrange=[ydownvel,yupvel],ystyle=1

if (n_elements(color) ge n1) then begin
    loadct,4
    for j=0,n1-1 do begin
        oplot,velarr[j,*],psym=symbol[j],color=color[j]
        xyouts,0.1,0.9-(j*0.03),datalist[j],color=color[j],/normal
    endfor
endif else begin
    for j=0,n1-1 do begin
        oplot,velarr[j,*],psym=symbol[j],color=60+(j*10)
        xyouts,0.1,0.9-(j*0.03),datalist[j],color=60+(j*10),/normal
    endfor
endelse

pause
window,0,retain=2,xsize=900
device,decomposed=0

loadct,0
plot,sigarr[0,*],/ynozero,/nodata,title='VELOCITY DISPERSION',$
  yrange=[ydownsig,yupsig],ystyle=1

if (n_elements(color) ge n1) then begin
    loadct,4
    for j=0,n1-1 do begin
        oplot,sigarr[j,*],psym=symbol[j],color=color[j]
        xyouts,0.1,0.9-(j*0.03),datalist[j],color=color[j],/normal
    endfor
endif else begin
    for j=0,n1-1 do begin
        oplot,sigarr[j,*],psym=symbol[j],color=60+(j*10)
        xyouts,0.1,0.9-(j*0.03),datalist[j],color=60+(j*10),/normal
    endfor
endelse
stop
end
