pro binmaker2

readcol,'F2_D1.coord',f='i,f,f,f',fib,rad,ra,dec
n0= n_elements(ra)

readcol,'bin.range',f='f,f',r1,r2
n1 = n_elements(r1)
;colors = intarr(n1)
;for j=0,n1-1 do colors[j] = (j * floor(165.0/n1))
;print,colors & pause
colors = [20,60,110,180,150,255]

window,3,retain=2
device,decomposed=0
loadct,0
plot,ra,dec,/nodata,/ynozero
for j=0,245 do xyouts,ra[j],dec[j],strn(fib[j])
loadct,4

irad = bsort(rad)
srad = rad[irad]
sra = ra[irad]
sdec = dec[irad]
sfib = fib[irad]

for j=0,n1-1 do begin
    ibin = where(srad gt r1[j] and srad le r2[j])
    print,r1[j],' ',r2[j],' ',n_elements(ibin)
    radp = srad[ibin]
    n2 = n_elements(radp)
    for k=0,n2-1 do xyouts,sra[ibin[k]],sdec[ibin[k]],strn(sfib[ibin[k]]),$
      color=colors[j]
    binname = ''
    print,'Enter a bin name...'
    read,binname
    free_lun,5
    openw,5,binname
    for k=0,n2-1 do begin
        a1 = strn(sfib[ibin[k]])+'_d1'
        a2 = strn(sfib[ibin[k]])+'_d2'
        a3 = strn(sfib[ibin[k]])+'_d3'
        printf,5,a1
        printf,5,a2
        printf,5,a3
    endfor
    free_lun,5
;    for k=0,n2-1 do plots,sra[ibin[k]],sdec[ibin[k]],psym=1,$
;      symsize=1.5,thick=2,color=colors[j]
    pause
endfor

;for j=0,245 do begin
;    openw,5,'M32d3_'+strn(j+1)		
;    printf,5,strn(j+1)+'_d3'
;    free_lun,5
;endfor

stop
end
