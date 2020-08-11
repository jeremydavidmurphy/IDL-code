; This routine uses the pfitlov.out file and plots a velocity,
; velocity dispersion, h3 and h4 profile. THIS DOES NOT USE THE
; CONTOUR FUNCTION but rather returns dots. It's been folded into
; Pfield.pro.

pro justplot, file

readcol,file,format='x,x,f,f,f,f,f,f',ra,dec,vel,disp,h3,h4,silent=1

nnn = n_elements(ra)
rup = max(vel)
rdown = min(vel)

set_plot,'ps'
device,file='M87vel_dots.ps',/color
loadct,0
levels = (rup-rdown)/255.*findgen(255)+rdown
loadct,0

spotcolor = intarr(nnn)
si = bsort(vel)
sv = vel[si]
sra = ra[si]
sdec = dec[si]

slope = 255./(max(sv)-min(sv))
for j=0,nnn-1 do spotcolor[j] = slope * (sv[j]-min(sv))

loadct,0
plot,sra,sdec,psym=3,/nodata,/isotropic,position=[0.33,0.10,0.85,0.88],$
  xthick=3,ythick=3,xrange=[min(ra)-10,max(ra)+10],/xstyle,xtitle='arcsec',$
  ytitle='arcsec',charthick=2,yrange=[min(dec)-10,max(dec)+10],/ystyle
loadct,33
for j=0,nnn-1 do plots,sra[j],sdec[j],psym=sym(1),color=spotcolor[j],$
  symsize=0.7

ncolors = n_elements(spotcolor)
loc = [0.27,0.915,0.63,0.945]
bar = bindgen(256) # replicate(1b,10)
xsize = (loc[2] - loc[0]) * !d.x_vsize
ysize = (loc[3] - loc[1]) * !d.y_vsize
xstart = loc[0] * !d.x_vsize
ystart = loc[1] * !d.y_vsize
tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
loadct,0
plots, [loc[0], loc[0], loc[2], loc[2], loc[0]], [loc[1], loc[3], loc[3], loc[1], loc[1]],/normal

middle = (max(vel)-min(vel))/2
vup = middle
vdown = -middle
xyouts,0.27,0.89,strn(vdown),/normal,charthick=2,charsize=0.8
xyouts,0.565,0.89,strn(vup),/normal,charthick=2,charsize=0.8
xyouts,0.42,0.89,'km/sec',/normal,charthick=2,charsize=0.8
xyouts,0.40,0.955,'M87 Velocity',/normal,charthick=2,charsize=1.0
device,/close_file

;******************************************************************
jump:
rup = max(disp)
rdown = min(disp)
set_plot,'ps'
device,file='M87disp_dots.ps'
loadct,0

spotcolor = intarr(nnn)
si = bsort(disp)
sv = disp[si]
sra = ra[si]
sdec = dec[si]

slope = 230./(max(sv)-min(sv))
for j=0,nnn-1 do spotcolor[j] = slope * (sv[j]-min(sv))

loadct,0
plot,sra,sdec,psym=3,/nodata,/isotropic,position=[0.33,0.10,0.85,0.88],$
  xthick=3,ythick=3,xrange=[min(ra)-10,max(ra)+10],/xstyle,xtitle='arcsec',$
  ytitle='arcsec',charthick=2,yrange=[min(dec)-10,max(dec)+10],/ystyle
loadct,33
for j=0,nnn-1 do plots,sra[j],sdec[j],psym=sym(1),color=spotcolor[j],$
  symsize=0.7

loc = [0.27,0.915,0.63,0.945]
bar = bindgen(240) # replicate(1b,10)
xsize = (loc[2] - loc[0]) * !d.x_vsize
ysize = (loc[3] - loc[1]) * !d.y_vsize
xstart = loc[0] * !d.x_vsize
ystart = loc[1] * !d.y_vsize
tv, bar, xstart, ystart, xsize=xsize, ysize=ysize
loadct,0
plots, [loc[0], loc[0], loc[2], loc[2], loc[0]], [loc[1], loc[3], loc[3], loc[1], loc[1]],/normal

;middle = (max(disp)-min(disp))/2
;vup = middle
;vdown = -middle
dup = max(disp)
ddown = min(disp)
xyouts,0.27,0.89,strn(ddown),/normal,charthick=2,charsize=0.8
xyouts,0.565,0.89,strn(dup),/normal,charthick=2,charsize=0.8
xyouts,0.42,0.89,'km/sec',/normal,charthick=2,charsize=0.8
xyouts,0.39,0.955,'M87 Dispersion',/normal,charthick=2,charsize=1.0
device,/close_file

end
