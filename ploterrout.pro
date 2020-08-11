PRO ploterrout, list, radiusfile, scale=scale

; This routine just plots the moments, with uncertainties, from the
; err.out files. This is the last step run in the extraction of the
; kinematics and expresses the best modeling moments and their
; uncertainties. What's plotted is the biweight and biweight
; uncertainty of each spectral region and their 100 MC simulations

; LIST: A list of all the bin####_err.out files you want to plot.

; RADIUSFILE: A list of 2 columns of the following format:
; 01     -7.002
; 05     -6.417
; 02    -11.912
; 05     -8.919
; 01    -14.254
; 04    -12.814
; Where the first column is the angular bin and the second is the
; radial distance, in arcseconds.
 
; SCALE: This is the scale of the galaxy (parsecs/arcseconds). If you
; include this, the 2nd x-axis is given in units of kpc.

; The output is a ps figure showing the four moments.

yv1 = -120
yv2 = 120
yd1 = 225
yd2 = 350
yh31 = -0.1
yh32 = 0.1
yh41 = -0.1
yh42 = 0.1


readcol,list,silent=1,f='a',files
readcol,radiusfile,silent=1,f='i,f',bins,radius
i1 = where(bins eq 01)
i2 = where(bins eq 02)
i3 = where(bins eq 03)
i4 = where(bins eq 04)
i5 = where(bins eq 05)
colors = [60,110,180,150]

n0 = n_elements(files)
if (n0 ne n_elements(radius)) then begin
    print,'Your lists are not the same length!'
    stop
endif

moments = fltarr(4,n0,3) ;moments, bins, -uncertainty-to-+uncertainty

for j=0,n0-1 do begin
    readcol,files[j],silent=1,f='x,f,f,x,f,f,x,f,f,x,f,f',a,b,c,d,e,f,g,h
    moments[0,j,0] = a
    moments[1,j,0] = c
    moments[2,j,0] = e
    moments[3,j,0] = g
    moments[0,j,1] = a-b
    moments[1,j,1] = c-d
    moments[2,j,1] = e-f
    moments[3,j,1] = g-h
    moments[0,j,2] = a+b
    moments[1,j,2] = c+d
    moments[2,j,2] = e+f
    moments[3,j,2] = g+h
endfor

x1 = min(radius)-25.0
x2 = max(radius)+25.0
set_plot,'ps'
device,file='moments_w_errors.ps',/color
!p.region = [0.1,0.1,0.9,0.9]
!p.multi = [0,1,4,0,1]
!y.omargin = [2,4]

;VELOCITY
if (n_elements(scale) ne 0) then begin
plot,radius[i1],moments[0,i1,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
  xrange=[x1,x2],xstyle=9,yrange=[yv1,yv2],/ys,pos=[0.15,0.675,0.85,0.85],$
  xtickformat='(A1)',xticklen=0.07,yminor=2,/nodata,symsize=1.5
xyouts,0.105,0.7,'Velocity (km/s)',/normal,orientation=90,charsize=0.6,$
  charthick=1.5
endif else begin
plot,radius[i1],moments[0,i1,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
  xrange=[x1,x2],xstyle=1,yrange=[yv1,yv2],/ys,pos=[0.15,0.675,0.85,0.85],$
  xtickformat='(A1)',xticklen=0.07,yminor=2,/nodata,symsize=1.5
xyouts,0.105,0.7,'Velocity (km/s)',/normal,orientation=90,charsize=0.6,$
  charthick=1.5
endelse
loadct,4,/silent
plots,radius[i2],moments[0,i2,0],psym=sym(9),color=60,thick=3,symsize=0.9
plots,radius[i3],moments[0,i3,0],psym=sym(9),color=110,thick=3,symsize=0.9
plots,radius[i4],moments[0,i4,0],psym=sym(9),color=180,thick=3,symsize=0.9
plots,radius[i5],moments[0,i5,0],psym=sym(9),color=150,thick=3,symsize=0.9
loadct,0,/silent
plots,radius[i1],moments[0,i1,0],psym=sym(4),thick=3,symsize=0.9
for j=0,n_elements(i1)-1 do oplot,[radius[i1[j]],radius[i1[j]]],$
  [moments[0,i1[j],1],moments[0,i1[j],2]],thick=2

if (n_elements(scale) ne 0) then begin ;the second xaxis is put in kpc)
    irad = bsort(radius)
    radkpc = (radius * scale) / 1000.0
    radkpc = radkpc[irad]
    x12 = min(radkpc)-((25.0*scale)/1000.0)
    x22 = max(radkpc)+((25.0*scale)/1000.0)
loadct,0,/silent
axis,xaxis=1,xrange=[x12,x22],/save,xtitle='Radius (kpc)',$
  charthick=2,xthick=2,charsize=1.3,xstyle=1,xticklen=0.07
endif

;VELOCITY DISPERSION
plot,radius[i1],moments[1,i1,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
  xrange=[x1,x2],/xs,yrange=[yd1,yd2],/ys,pos=[0.15,0.50,0.85,0.675],$
  xtickformat='(A1)',xticklen=0.07,yminor=2,/nodata,symsize=1.5
xyouts,0.105,0.51,'Dispersion (km/s)',/normal,orientation=90,charsize=0.6,$
  charthick=1.5
loadct,4,/silent
plots,radius[i2],moments[1,i2,0],psym=sym(9),color=60,thick=3,symsize=0.9
plots,radius[i3],moments[1,i3,0],psym=sym(9),color=110,thick=3,symsize=0.9
plots,radius[i4],moments[1,i4,0],psym=sym(9),color=180,thick=3,symsize=0.9
plots,radius[i5],moments[1,i5,0],psym=sym(9),color=150,thick=3,symsize=0.9
loadct,0,/silent
plots,radius[i1],moments[1,i1,0],psym=sym(4),thick=3,symsize=0.9
for j=0,n_elements(i1)-1 do oplot,[radius[i1[j]],radius[i1[j]]],$
  [moments[1,i1[j],1],moments[1,i1[j],2]],thick=2

;H3
plot,radius[i1],moments[2,i1,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
  xrange=[x1,x2],/xs,yrange=[yh31,yh32],/ys,pos=[0.15,0.325,0.85,0.50],$
  xtickformat='(A1)',xticklen=0.07,yminor=2,/nodata,symsize=1.5
xyouts,0.105,0.412,'H3',/normal,orientation=90,charsize=0.6,$
  charthick=1.5
loadct,4,/silent
plots,radius[i2],moments[2,i2,0],psym=sym(9),color=60,thick=3,symsize=0.9
plots,radius[i3],moments[2,i3,0],psym=sym(9),color=110,thick=3,symsize=0.9
plots,radius[i4],moments[2,i4,0],psym=sym(9),color=180,thick=3,symsize=0.9
plots,radius[i5],moments[2,i5,0],psym=sym(9),color=150,thick=3,symsize=0.9
loadct,0,/silent
plots,radius[i1],moments[2,i1,0],psym=sym(4),thick=3,symsize=0.9
for j=0,n_elements(i1)-1 do oplot,[radius[i1[j]],radius[i1[j]]],$
  [moments[2,i1[j],1],moments[2,i1[j],2]],thick=2
oplot,[-1000,1000],[0.0,0.0],linestyle=1

;H4
plot,radius[i1],moments[3,i1,0],xthick=2,ythick=2,charthick=2,psym=sym(6),$
  xrange=[x1,x2],/xs,yrange=[yh41,yh42],/ys,pos=[0.15,0.15,0.85,0.325],$
  xticklen=0.07,yminor=2,/nodata,symsize=1.5,xtitle='Radius (arcsec)',$
  charsize=1.0,xcharsize=1.3
xyouts,0.105,0.237,'H4',/normal,orientation=90,charsize=0.6,$
  charthick=1.5
loadct,4,/silent
plots,radius[i2],moments[3,i2,0],psym=sym(9),color=60,thick=3,symsize=0.9
plots,radius[i3],moments[3,i3,0],psym=sym(9),color=110,thick=3,symsize=0.9
plots,radius[i4],moments[3,i4,0],psym=sym(9),color=180,thick=3,symsize=0.9
plots,radius[i5],moments[3,i5,0],psym=sym(9),color=150,thick=3,symsize=0.9
loadct,0,/silent
plots,radius[i1],moments[3,i1,0],psym=sym(4),thick=3,symsize=0.9
for j=0,n_elements(i1)-1 do oplot,[radius[i1[j]],radius[i1[j]]],$
  [moments[3,i1[j],1],moments[3,i1[j],2]],thick=2
oplot,[-1000,1000],[0.0,0.0],linestyle=1

device,/close_file

stop
END
