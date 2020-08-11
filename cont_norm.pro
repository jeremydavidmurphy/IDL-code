PRO cont_norm, list

box = 150

data = readfits('N4472d.fits')
dsmtd = smooth(data,box,/edge_truncate,/nan)
sdata = data / dsmtd
 
readcol,list,f='a',files
n0 = n_elements(files)
outout = fltarr(2048,n0)
window,3,retain=2,xsize=600,ysize=300
window,2,retain=2,xsize=600,ysize=300
device,decomposed=0

for j=0,n0-1 do begin
    spec = readfits(files[j])
    wave = spec[*,0]
    spec = spec[*,1]
    test = median(spec) * 2.0
    i = where(spec gt test,count)
    if (count gt 0) then begin
        sspec = spec
        sspec[i] = !Values.F_NAN
    endif else sspec = spec
    smtd = smooth(sspec,box,/edge_truncate,/nan)
    
    out = spec / smtd
    
    wset,2
    loadct,0
    plot,wave,spec,yrange=[0,100]
    loadct,4
    oplot,wave,smtd,color=255
    wset,3
;    if (j eq 0) then begin
    plot,wave,out,yrange=[0,2]
;    endif
;    oplot,wave,out,color=110
    wait,2
    outout[*,j] = out
endfor


var = fltarr(2048)
for j=0,2047 do var[j] = variance(outout[j,*])
stddev = sqrt(var)
wset,2
loadct,0
plot,wave,var,yrange=[0,0.01]
pause
plot,wave,stddev

stop
END
