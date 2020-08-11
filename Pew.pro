PRO Pew, list

readcol,list,f='a',ewlist
n0 = n_elements(ewlist)

readcol,'bin.radius',f='x,f',radius
if (n_elements(radius) ne n0) then stop

readcol,ewlist[0],f='x,f,x,a',test,names
n1 = n_elements(test)

ewarr = fltarr(n0,n1)

;readcol,'eqw.dat',f='x,x,x,x,x,x,x,a',spec
;readcol,'bin.index',f='i',index

for j=0,n0-1 do begin
    readcol,ewlist[j],f='x,f,x,x',ew
    ewarr[j,*] = ew
endfor

;colors = [60,110,180,150]

window,0,retain=2
device,decomposed=0
for j=0,n1-1 do begin
    loadct,0
    plot,radius,ewarr[*,j],psym=2,/ynozero,title=names[j],charsize=1.2,$
      xtitle='Radius (arcsec)',ytitle='EW (A)'
;    i = where(index eq 1,count)
;    if (count gt 0) then oplot,radius[i],ewarr[i,j],psym=2
;    loadct,4
;    for k=1,4 do begin
;        i = where(index eq k+1,count)
;        if (count gt 0) then oplot,radius[i],ewarr[i,j],psym=2,color=colors[k-1]
;    endfor
    pause
;    if (j eq n1-1) then wdelete
endfor
stop
set_plot,'ps'
for j=0,n1-1 do begin
    device,file=spec[j]+'.ps',/color
    loadct,0
    plot,radius,ewarr[*,j],psym=2,/ynozero,title=spec[j],$
      charsize=1.0,xthick=3,ythick=3,charthick=3,thick=3,$
      xtitle='Radius (arcsec)',ytitle='EW (A)',/nodata
    i = where(index eq 1,count)
    if (count gt 0) then oplot,radius[i],ewarr[i,j],psym=2,thick=3
    loadct,4
    for k=1,4 do begin
        i = where(index eq k+1,count)
        if (count gt 0) then oplot,radius[i],ewarr[i,j],psym=2,color=colors[k-1],thick=3
    endfor
    device,/close_file
endfor
set_plot,'x'

stop
END
