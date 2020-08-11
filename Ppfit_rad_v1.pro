PRO Ppfit_rad, file, GALNAME=galname

; This routine is used to plot the pfitlov.out values AGAINST
; RADIUS. It requires the file, 'bin.radius' to exist. This file is of
; the form
;         bin####.coord   4.560
; where the second column are the radial coordinates. 

; FILE: The name of the pfitlov.out file

; GALNAME: If you don't enter this, you will be prompted to do so by
; the code. This is just for titles in the figures and naming of the
; output files.

; ERROR: Not working at present. Will allow a file of error values to
; be incorporated in the near future. This will then include error
; bars on all the moments.

if (n_elements(galname) eq 0) then begin
    galname = ''
    print,'Enter the name of the galaxy:'
    read,galname
endif

readcol,file,format='a,f,f,x,f,f,x,x,x',bin1,vel,disp,h3,h4
n0 = n_elements(bin1)
readcol,'bin.radius',format='a,f',bin2,rad
n1 = n_elements(bin2)

;if (n_elements(bin2) ne n0) then begin
;    print,'Your bin.radius and '+file+' are not the same length!'
;    print,'Fix this

bin1t = strarr(n0)
bin2t = strarr(n1)

for j=0,n0-1 do begin
    t = strsplit(bin1[j],'.',/extract)
    bin1t[j] = t[0]
endfor
for j=0,n1-1 do begin
    t = strsplit(bin2[j],'.',/extract)
    bin2t[j] = t[0]
endfor

i = where(bin1t eq bin2t)
binp = bin1t[i]
n2 = n_elements(binp)

radt = rad[i]

; Now the log, + and - is established
radlin = fltarr(n2)
radlog = fltarr(n2)
icolor = intarr(n2)

for j=0,n2-1 do begin
    t = strsplit(binp[j],'0123456789',/extract)
    tt = strsplit(binp[j],'n0',/extract)
    if (t[0] eq 'bin') then radlog[j] = alog10(radt[j])
    if (t[0] eq 'bin') then radlin[j] = radt[j]
    if (t[0] eq 'bnn') then radlog[j] = alog10(radt[j]) * (-1.0)
    if (t[0] eq 'bnn') then radlin[j] = radt[j] * (-1.0)
    if (tt[2] eq 1) then icolor[j] = 0
    if (tt[2] eq 2) then icolor[j] = 60
    if (tt[2] eq 3) then icolor[j] = 110
    if (tt[2] eq 4) then icolor[j] = 180
    if (tt[2] eq 5) then icolor[j] = 150
endfor

xup1 = max(radlin[i])+20.0
xdn1 = min(radlin[i])-20.0
xup2 = max(radlog[i])+0.2
xdn2 = min(radlog[i])-0.2

set_plot,'ps'
device,file=galname+'.v.ps',/color

!p.multi = [0,2,2,0,1]
!y.omargin = [2,4]

loadct,0
plot,radlin,vel[i],psym=1,title='Velocity for '+galname,ytitle='Velocity (km/sec)',$
  xtitle='Radius (arcsec)',charsize=0.8,xthick=3,ythick=3,thick=3,/ynozero,$
  charthick=3,xrange=[xdn1,xup1],xstyle=1,/nodata
loadct,4
for j=0,n2-1 do plots,radlin[j],vel[j],psym=1,color=icolor[j],thick=3

loadct,0
plot,radlog,vel[i],psym=1,title='Log Velocity for '+galname,ytitle='Velocity (km/sec)',$
  xtitle='Log Radius (arcsec)',charsize=0.8,xthick=3,ythick=3,thick=3,/ynozero,$
  charthick=3,xrange=[xdn2,xup2],xstyle=1,/nodata
loadct,4
for j=0,n2-1 do plots,radlog[j],vel[j],psym=1,color=icolor[j],thick=3

loadct,0
plot,radlin,disp[i],psym=1,title='Velocity Dispersion for '+galname,ytitle='Velocity Dispersion (km/sec)',$
  xtitle='Radius (arcsec)',charsize=0.8,xthick=3,ythick=3,thick=3,/ynozero,$
  charthick=3,xrange=[xdn1,xup1],xstyle=1,/nodata
loadct,4
for j=0,n2-1 do plots,radlin[j],disp[j],psym=1,color=icolor[j],thick=3

loadct,0
plot,radlog,disp[i],psym=1,title='Log Velocity Dispersion for '+galname,ytitle='Velocity Dispersion (km/sec)',$
  xtitle='Log Radius (arcsec)',charsize=0.8,xthick=3,ythick=3,thick=3,/ynozero,$
  charthick=3,xrange=[xdn2,xup2],xstyle=1,/nodata
loadct,4
for j=0,n2-1 do plots,radlog[j],disp[j],psym=1,color=icolor[j],thick=3

device,/close_file

device,file=galname+'.h.ps'

!p.multi = [0,2,2,0,1]
!y.omargin = [2,4]

loadct,0
plot,radlin,h3[i],psym=1,title='H3 for '+galname,ytitle='H3',$
  xtitle='Radius (arcsec)',charsize=0.8,xthick=3,ythick=3,thick=3,/ynozero,$
  charthick=3,xrange=[xdn1,xup1],xstyle=1,/nodata
loadct,4
for j=0,n2-1 do plots,radlin[j],h3[j],psym=1,color=icolor[j],thick=3

loadct,0
plot,radlog,h3[i],psym=1,title='Log H3 for '+galname,ytitle='H3',$
  xtitle='Log Radius (arcsec)',charsize=0.8,xthick=3,ythick=3,thick=3,/ynozero,$
  charthick=3,xrange=[xdn2,xup2],xstyle=1,/nodata
loadct,4
for j=0,n2-1 do plots,radlog[j],h3[j],psym=1,color=icolor[j],thick=3

loadct,0
plot,radlin,h4[i],psym=1,title='H4 for '+galname,ytitle='H4',$
  xtitle='Radius (arcsec)',charsize=0.8,xthick=3,ythick=3,thick=3,/ynozero,$
  charthick=3,xrange=[xdn1,xup1],xstyle=1,/nodata
loadct,4
for j=0,n2-1 do plots,radlin[j],h4[j],psym=1,color=icolor[j],thick=3

loadct,0
plot,radlog,h4[i],psym=1,title='Log H4 for '+galname,ytitle='H4',$
  xtitle='Log Radius (arcsec)',charsize=0.8,xthick=3,ythick=3,thick=3,/ynozero,$
  charthick=3,xrange=[xdn2,xup2],xstyle=1,/nodata
loadct,4
for j=0,n2-1 do plots,radlog[j],h4[j],psym=1,color=icolor[j],thick=3

device,/close_file
set_plot,'x'

stop
END
