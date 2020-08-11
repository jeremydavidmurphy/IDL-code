PRO Ppfitn2832,type


readcol,'bin.radius.N2832',f='a,f',bins,rad
readcol,'pfit.list',f='a',lists
n0 = n_elements(lists)

readcol,'bin.colors',f='i',silent=1,colors2


readcol,lists[0],f='a,f,f,x,f,f,x,x,x',names,v,d,h3,h4,silent=1
colors = fltarr(n0)
for j=0,n0-1 do colors[j] = j * (255.0/n0)
colors = floor(colors)
colors = colors + colors[1]

if (type eq 'screen') then begin
set_plot,'x'
window,0,retain=2,xsize=1200
device,decomposed=0

loadct,0,/silent
plot,rad,v,psym=1,/nodata,xrange=[-150,150],yrange=[-100,100],$
  /xs,/ys,title='Velocity'
loadct,27,/silent

for j=0,n0-1 do begin
    readcol,lists[j],f='x,f,f,x,f,f,x,x,x',v,d,h3,h4,silent=1
    oplot,rad,v,psym=1,color=colors[j]
    pause
endfor

loadct,0,/silent
plot,rad,v,psym=1,/nodata,xrange=[-150,150],yrange=[200,400],$
  /xs,/ys,title='Velocity Dispersion'
loadct,27,/silent

for j=0,n0-1 do begin
    readcol,lists[j],f='x,f,f,x,f,f,x,x,x',v,d,h3,h4,silent=1
    oplot,rad,d,psym=1,color=colors[j]
    pause
endfor

wdelete,0
endif

if (type eq 'ps') then begin
arr = fltarr(77,n0,4)
set_plot,'ps'
device,file='N2832_vel.ps',/color

loadct,0,/silent
plot,rad,v,psym=1,/nodata,xrange=[-150,150],yrange=[-100,100],$
  /xs,/ys,title='Velocity',xthick=2,ythick=2,charthick=2,$
  xtitle='Radius (arcsec)',ytitle='Rotational Velocity (km/sec)'
loadct,27,/silent

for j=0,n0-1 do begin
    readcol,lists[j],f='x,f,f,x,f,f,x,x,x',v,d,h3,h4,silent=1
    arr[*,j,0] = v
    arr[*,j,1] = d
    arr[*,j,2] = h3
    arr[*,j,3] = h4
    oplot,rad,v,psym=sym(9),color=colors[j],thick=2
endfor

medvel = median(arr[*,*,0],dim=2,/even)
loadct,0,/silent
i = where(colors2 eq 0)
oplot,rad[i],medvel[i],psym=sym(1)
loadct,4,/silent
i = where(colors2 eq 1)
oplot,rad[i],medvel[i],psym=sym(1),color=60
i = where(colors2 eq 2)
oplot,rad[i],medvel[i],psym=sym(1),color=110
i = where(colors2 eq 3)
oplot,rad[i],medvel[i],psym=sym(1),color=250
i = where(colors2 eq 4)
oplot,rad[i],medvel[i],psym=sym(1),color=180
i = where(colors2 eq 5)
oplot,rad[i],medvel[i],psym=sym(1),color=150

device,/close_file

device,file='N2832_disp.ps',/color

loadct,0,/silent
plot,rad,h3,psym=1,/nodata,xrange=[-150,150],yrange=[200,400],$
  /xs,/ys,title='Velocity Dispersion',xthick=2,ythick=2,charthick=2,$
  xtitle='Radius (arcsec)',ytitle='Velocity Dispersion (km/sec)'
loadct,27,/silent

for j=0,n0-1 do oplot,rad,arr[*,j,1],psym=sym(9),color=colors[j],thick=2

medvel = median(arr[*,*,1],dim=2,/even)
loadct,0,/silent
oplot,rad,medvel,psym=sym(1)
device,/close_file

device,file='N2832_h3.ps',/color

loadct,0,/silent
plot,rad,h3,psym=1,/nodata,xrange=[-150,150],yrange=[-0.5,0.5],$
  /xs,/ys,title='Gauss-Hermite Coefficient 3',xthick=2,ythick=2,charthick=2,$
  xtitle='Radius (arcsec)',ytitle='H3'
loadct,27,/silent

for j=0,n0-1 do oplot,rad,arr[*,j,2],psym=sym(9),color=colors[j],thick=2

medvel = median(arr[*,*,2],dim=2,/even)
loadct,0,/silent
oplot,rad,medvel,psym=sym(1)
device,/close_file

device,file='N2832_h4.ps',/color

loadct,0,/silent
plot,rad,d,psym=1,/nodata,xrange=[-150,150],yrange=[-0.5,0.5],$
  /xs,/ys,title='Gauss-Hermite Coefficient 4',xthick=2,ythick=2,charthick=2,$
  xtitle='Radius (arcsec)',ytitle='H4'
loadct,27,/silent

for j=0,n0-1 do oplot,rad,arr[*,j,3],psym=sym(9),color=colors[j],thick=2

medvel = median(arr[*,*,3],dim=2,/even)
loadct,0,/silent
oplot,rad,medvel,psym=sym(1)
device,/close_file

set_plot,'x'
endif

stop
END
