PRO PpfitMED, list, outname=outname

; This routine is a modification of pfitave.pro and now works from a
; list. This list contains pfitlov.out files. The routine will read
; these in, in turn, and return median values for the 4 moments, and
; SD estimates for each.

; This runs parallel with Ppfit_rad.pro. It also reads in the
; bin.radius file 

readcol,list,f='a',silent=1,files
n0 = n_elements(files)

readcol,files[0],format='a,f,f,x,f,f,x,x,x',binname,v,d,h3,h4,silent=1
n1 = n_elements(v)

data = fltarr(4,n1,n0)

for j=0,n0-1 do begin
    readcol,files[j],format='a,f,f,x,f,f,x,x,x',binname,v,d,h3,h4
    data[0,*,j] = v
    data[1,*,j] = d
    data[2,*,j] = h3
    data[3,*,j] = h4
endfor

readcol,'bin.radius',f='x,f',radius

final = fltarr(4,n1,2)
for j=0,n1-1 do begin
    final[*,j,0] = median(data[*,j,*],dim=3)
    for k=0,3 do final[k,j,1] = stddev(data[k,j,*])
endfor

; Now, a cheap velocity correction is made
ismall = where(abs(radius) le 10.0)
voff = median(final[0,ismall,0])
final[0,*,0] = final[0,*,0] - voff

for j=0,n1-1 do begin
    trim = strsplit(binname[j],'_',/extract)
    binname[j] = trim[0]
endfor

binp = strarr(n1)

for j=0,n1-1 do begin
    t = strsplit(binname[j],'c_.',/extract)
    binp[j] = t[0]
endfor

i1 = 0
i2 = 0
i3 = 0
i4 = 0
i5 = 0

for j=0,n1-1 do begin
    t = strsplit(binp[j],'0123456789',/extract)
    tt = strsplit(binp[j],'n0',/extract)
    if (tt[2] eq 1) then i1 = [i1,j]
    if (tt[2] eq 2) then i2 = [i2,j]
    if (tt[2] eq 3) then i3 = [i3,j]
    if (tt[2] eq 4) then i4 = [i4,j]
    if (tt[2] eq 5) then i5 = [i5,j]
endfor
i1 = i1[1:*]
i2 = i2[1:*]
i3 = i3[1:*]
i4 = i4[1:*]
i5 = i5[1:*]
icolor = [60,110,180,150]

eup1 = final[0,*,0] + final[0,*,1]
edn1 = final[0,*,0] - final[0,*,1]
eup2 = final[1,*,0] + final[1,*,1]
edn2 = final[1,*,0] - final[1,*,1]

window,0,retain=2
device,decomposed=0
loadct,0
plot,radius,final[0,*,0],psym=2
errplot,radius[i1],edn1[i1],eup1[i1]
loadct,4
oplot,radius[i5],final[0,i5,0],color=icolor[3],psym=2
oplot,radius[i4],final[0,i4,0],color=icolor[2],psym=2
oplot,radius[i3],final[0,i3,0],color=icolor[1],psym=2
oplot,radius[i2],final[0,i2,0],color=icolor[0],psym=2
errplot,radius[i5],edn1[i5],eup1[i5],color=icolor[3]
errplot,radius[i4],edn1[i4],eup1[i4],color=icolor[2]
errplot,radius[i3],edn1[i3],eup1[i3],color=icolor[1]
errplot,radius[i2],edn1[i2],eup1[i2],color=icolor[0]
loadct,0
oplot,radius[i1],final[0,i1,0],psym=2
errplot,radius[i1],edn1[i1],eup1[i1]

if (n_elements(outname) eq 0) then begin
    outname = ''
    print,'Enter the name of the output file:'
    read,outname
    out1 = 'Vel_'+outname+'.ps'
    out2 = 'Disp_'+outname+'.ps'
endif else begin
    out1 = 'Vel_'+outname+'.ps'
    out2 = 'Disp_'+outname+'.ps'
endelse

set_plot,'ps'
device,file=out1,/color

loadct,0
plot,radius,final[0,*,0],psym=4,title='N4472 Velocity Profile- '+outname,$
  xtitle='Radius (arcsec)',ytitle='Velocity (km/sec)',xthick=2,ythick=2,$
  charthick=2,yrange=[-150,150],/ystyle,thick=2
errplot,radius[i1],edn1[i1],eup1[i1],thick=2
loadct,4
oplot,radius[i5],final[0,i5,0],color=icolor[3],psym=4,thick=2
oplot,radius[i4],final[0,i4,0],color=icolor[2],psym=4,thick=2
oplot,radius[i3],final[0,i3,0],color=icolor[1],psym=4,thick=2
oplot,radius[i2],final[0,i2,0],color=icolor[0],psym=4,thick=2
errplot,radius[i5],edn1[i5],eup1[i5],color=icolor[3],thick=2
errplot,radius[i4],edn1[i4],eup1[i4],color=icolor[2],thick=2
errplot,radius[i3],edn1[i3],eup1[i3],color=icolor[1],thick=2
errplot,radius[i2],edn1[i2],eup1[i2],color=icolor[0],thick=2
loadct,0
oplot,radius[i1],final[0,i1,0],psym=4,thick=2
errplot,radius[i1],edn1[i1],eup1[i1],thick=2
;xyouts,0.2,0.2,strn(voff),/normal

device,/close_file

pause
set_plot,'x'
window,0,retain=2
device,decomposed=0
loadct,0
plot,radius,final[1,*,0],psym=4,yrange=[200,600]
errplot,radius[i1],edn2[i1],eup2[i1]
loadct,4
oplot,radius[i5],final[1,i5,0],color=icolor[3],psym=4
oplot,radius[i4],final[1,i4,0],color=icolor[2],psym=4
oplot,radius[i3],final[1,i3,0],color=icolor[1],psym=4
oplot,radius[i2],final[1,i2,0],color=icolor[0],psym=4
errplot,radius[i5],edn2[i5],eup2[i5],color=icolor[3]
errplot,radius[i4],edn2[i4],eup2[i4],color=icolor[2]
errplot,radius[i3],edn2[i3],eup2[i3],color=icolor[1]
errplot,radius[i2],edn2[i2],eup2[i2],color=icolor[0]
loadct,0
oplot,radius[i1],final[1,i1,0],psym=4
errplot,radius[i1],edn2[i1],eup2[i1]

set_plot,'ps'
device,file=out2,/color

loadct,0
plot,radius,final[1,*,0],psym=4,title='N4472 Velocity Dispersion Profile- '+outname,$
  xtitle='Radius (arcsec)',ytitle='Velocity Dispersion (km/sec)',xthick=2,ythick=2,$
  charthick=2,yrange=[200,600],/ystyle,thick=2
errplot,radius[i1],edn2[i1],eup2[i1],thick=2
loadct,4
oplot,radius[i5],final[1,i5,0],color=icolor[3],psym=4,thick=2
oplot,radius[i4],final[1,i4,0],color=icolor[2],psym=4,thick=2
oplot,radius[i3],final[1,i3,0],color=icolor[1],psym=4,thick=2
oplot,radius[i2],final[1,i2,0],color=icolor[0],psym=4,thick=2
errplot,radius[i5],edn2[i5],eup2[i5],color=icolor[3],thick=2
errplot,radius[i4],edn2[i4],eup2[i4],color=icolor[2],thick=2
errplot,radius[i3],edn2[i3],eup2[i3],color=icolor[1],thick=2
errplot,radius[i2],edn2[i2],eup2[i2],color=icolor[0],thick=2
loadct,0
oplot,radius[i1],final[1,i1,0],psym=4,thick=2
errplot,radius[i1],edn2[i1],eup2[i1],thick=2
;youts,0.2,0.2,strn(voff),/normal

device,/close_file
set_plot,'x'

stop
END

