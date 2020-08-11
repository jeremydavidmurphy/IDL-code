PRO Pdensity

; This routine reads in dL.dat, dMhalo.dat, and dM.dat and then
; generates a density plot for a given galaxy. It also needs
; bindemo_r.out and bindemo_v.out in the calling directory.

; This code returns the density
ML = 7.0
rmax = 2250.0 ;The radial scaling factor
scale = 86.5 ;The number of pc per arcsec
xr = [0.2,rmax] ;the x-plotting range
yr = [10e-8,1000]
;*********************************************************************

readcol,'dL.dat',f='i,i,d',silent=1,rd,ad,dml
readcol,'dM.dat',f='x,x,d',silent=1,dmm
readcol,'dMhalo.dat',f='x,x,d',silent=1,dmh
readcol,'bindemo_r.out',silent=1,skipline=5,f='i,i,f,f,f,x,x',rr,ar,rl,rm,ru
readcol,'bindemo_v.out',silent=1,skipline=4,f='i,i,f,f,f',rv,av,al,am,au

readcol,'density.out',silent=1,f='x,x,x,x,x,d',dl
nr = n_elements(rr)
na = n_elements(rv)

dml = dml * ML

density = dblarr(na,nr,4);angular,radial,D3:density,L_density,DM_density,TOTAL_density
vol = dblarr(nr)

for j=0,nr-1 do begin
    r2 = ru[j] * rmax * scale ;the upper range in radius
    r1 = rl[j] * rmax * scale ;the lower range in radius
    vol[j] = double((!pi/15.0) * (r2^3 - r1^3))
    ir = where(rd eq rr[j])
    ml = dml[ir]
    mm = dmm[ir]
    mh = dmh[ir]
    for k=0,na-1 do density[k,j,0] = ml[k] / vol[j]
    for k=0,na-1 do density[k,j,1] = mh[k] / vol[j]
    for k=0,na-1 do density[k,j,2] = mm[k] / vol[j]
    density[*,j,3] = dl[ir]
;    if j eq 0 then stop
endfor

for k=0,na-1 do density[k,*,0] = smooth(density[k,*,0],3)
for k=0,na-1 do density[k,*,1] = smooth(density[k,*,1],3)
for k=0,na-1 do density[k,*,2] = smooth(density[k,*,2],3)

radius = rm * rmax
rkpc0 = (xr[0] * scale) / 1000.0
rkpc1 = (xr[1] * scale) / 1000.0

colors = intarr(na)
for j=0,na-1 do colors[j] = 60 + (j * floor(165.0/na))

window,0,retain=2
device,decomposed=0
plot,radius,density[0,*,0],/xlog,/ylog,xrange=[xr[0],xr[1]],/xstyle,$
  yrange=[yr[0],yr[1]]

loadct,4,/silent
for k=0,na-1 do begin
    oplot,radius,density[k,*,0],color=colors[k]
    oplot,radius,density[k,*,1],color=colors[k]
    oplot,radius,density[k,*,2],color=colors[k]
;    oplot,radius,density[k,*,3],color=colors[k]
endfor

set_plot,'ps'
device,file='density.ps',/color
loadct,0,/silent
plot,radius,density[0,*,0],/xlog,/ylog,xrange=[xr[0],xr[1]],xstyle=9,$
  yrange=[yr[0],yr[1]],xtitle='Radius (arcsec)',xthick=4,$
  ythick=4,charthick=4,ytitle='Density (M!dsun!n/pc!e3!n)',$
  position=[0.13,0.13,0.93,0.89],/ys,charsize=1.2

loadct,4,/silent
for k=0,na-1 do begin
    oplot,radius,density[k,*,0],color=colors[k]
    oplot,radius,density[k,*,2],color=colors[k]
endfor
loadct,0,/silent
oplot,radius,density[0,*,1],thick=5
axis,xaxis=1,xrange=[rkpc0,rkpc1],/save,xtitle='R (kpc)',$
  charthick=4,/xstyle,xthick=4,charsize=1.2

device,/close_file
set_plot,'x'

stop
END
