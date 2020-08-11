PRO blue_excess, list1, list2

readcol,list1,f='a',files1
n0 = n_elements(files1)
readcol,list2,f='a',files2

readcol,'bin.radius',f='x,f',radius
test = readfits(files1[0])
n1 = n_elements(test)

data = fltarr(n1,n0)
s2n = fltarr(n1,n0)
values = fltarr(7,n0)
ratio1 = fltarr(n0)
ratio2 = fltarr(n0)

wave = findgen(n1)*0.375 + 3550.0
for j=0,n0-1 do begin
    print,'Working on frame '+files1[j]+'...'
    data[*,j] = readfits(files1[j],/silent)
    s2n[*,j] = readfits(files2[j],/silent)
    values[0,j] = total(data[100:800,j])
    values[1,j] = median(data[100:800,j])
    values[2,j] = total(data[3000:4000,j])
    values[3,j] = median(data[3000:4000,j]);+75.0
    values[4,j] = median(s2n[100:800,j])
    values[5,j] = median(s2n[3000:4000,j])
    values[6,j] = median(s2n[*,j])
    ratio1[j] = values[1,j] / values[3,j]
    ratio2[j] = values[4,j] / values[5,j]
endfor

base = min(ratio1)
;base = ratio1[37]
for j=0,n0-1 do ratio2[j] = values[1,j] - values[3,j]*base

;imin = where(values[1,*] eq min(values[1,*]))
;v = ratio1[imin]
nratio1 = ratio1 / min(ratio1)

set_plot,'ps'

for j=0,n0-1 do begin
    device,file=files1[j]+'.ps'
    plot,wave,data[*,j],position=[0.12,0.55,0.94,0.93],ytitle='Data',$
      title='Spectra for '+files1[j],xthick=2,ythick=2,charthick=2,thick=2,$
      xrange=[wave[100],5500],/xstyle
    plot,wave,s2n[*,j],position=[0.12,0.1,0.94,0.47],/noerase,ytitle='S/N',$
      xthick=2,ythick=2,charthick=2,thick=2,xrange=[wave[100],5500],/xstyle
    xyouts,0.15,0.85,'Blue/Red ratio: '+strn(nratio1[j]),/normal,charthick=2,$
      charsize=1.3
    device,/close_file
endfor

set_plot,'x'
window,0,retain=2
device,decomposed=0
loadct,0
plot,radius,values[3,*]/min(values[3,*]),psym=1,/ylog,symsize=2,/nodata,$
yrange=[0.001,2],ystyle=9,position=[0.1,0.15,0.9,0.9],xrange=[-600,600]
loadct,4
oplot,radius,values[3,*]/max(values[3,*]),psym=2,symsize=2,color=150
oplot,radius,values[1,*]/max(values[1,*]),psym=4,symsize=2,color=60
;oplot,radius,values[3,*]/values[3,imin],psym=2,symsize=2,color=150
;oplot,radius,values[1,*]/values[1,imin],psym=4,symsize=2,color=60
;loadct,0
;oplot,radius,values[6,*]/max(values[6,*]),psym=5,symsize=2
;oplot,radius,values[6,*]/values[6,imin],psym=5,symsize=2
;loadct,4
xyouts,0.1,0.9,'Diamonds: 3590 - 3850 A flux',color=60,/normal
xyouts,0.1,0.85,'Stars: 4675 - 5050 A flux',color=150,/normal
xyouts,0.1,0.8,'Plus: Ratio of two values',color=110,/normal
axis,yaxis=1,yrange=[1,10],ytitle='Blue/Red Ratio (normalized to 1)',/save,color=110,/ystyle,xthick=2
;axis,yaxis=1,yrange=[min(nratio1),max(nratio1)],ytitle='Blue/Red Ratio (normalized to 1)',/save,color=110
oplot,radius,nratio1,psym=1,symsize=2,color=110

set_plot,'ps'
device,file='blue_excess.ps',/color
loadct,0
plot,radius,values[3,*]/min(values[3,*]),psym=1,/ylog,symsize=2,/nodata,$
  yrange=[0.001,2],ystyle=9,position=[0.13,0.15,0.9,0.9],xrange=[-400,400],$
  xthick=3,ythick=3,charthick=3,xtitle='Radius (arcsec)',$
  ytitle='Normalized Flux for Red and Blue Spectral Regions',$
  title='Data for N4472 and 1 pointing for M87'
loadct,4
oplot,radius,values[3,*]/max(values[3,*]),psym=2,symsize=2,color=150,thick=3
oplot,radius,values[1,*]/max(values[1,*]),psym=4,symsize=2,color=60,thick=3
xyouts,0.15,0.79,'Diamonds: 3590-3850A flux',color=60,/normal,charthick=3
xyouts,0.15,0.75,'Stars: 4675-5050A flux',color=150,/normal,charthick=3
xyouts,0.15,0.71,'Plus: Ratio of two values',color=110,/normal,charthick=3
oplot,[-76,-76],[0.001,2],thick=2
oplot,[76,76],[0.001,2],thick=2
loadct,4
xyouts,0.23,0.24,'M87f data',orientation=90,/normal,charthick=3,charsize=1.5
axis,yaxis=1,yrange=[1,5],ytitle='Blue/Red Ratio (normalized to 1)',/save,$
  /ystyle,ythick=3,charthick=3,/xstyle,color=110
oplot,radius,nratio1,psym=1,symsize=2,thick=5,color=110
device,/close

set_plot,'x'
loadct,0
plot,radius,values[3,*]/values[1,*],psym=2
oplot,[-76,-76],[0,1000]
oplot,[76,76],[0,1000]
pause
wdelete,0
stop
END
