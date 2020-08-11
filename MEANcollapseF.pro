function MEANcollapseF, data, maskname

; This routine runs a mean collapse on the data. This does a better
; job with the -666 flags, yet doesn't incldue the weighting.

aperture = 5.0
step = (aperture-1)/2.0

n1 = n_elements(data[*,0])

readcol,maskname,format='i,x',y
n3 = n_elements(y)
indexarr = y - 1.0

dataout = fltarr(n1,n3)

for jj=0,n3-1 do begin ;a loop through each fiber
    index = uint(indexarr[jj])
    rowD = data[*,index-step:index+step]
    for kk=0,n1-1 do begin ;a loop through wavelength
        top = rowD[kk,*]
        igood = where(top ne -666,count)
        if (count eq 5) then dataout[kk,jj] = mean(top)
        if (count gt 1) and (count le 4) then dataout[kk,jj] = mean(top[igood])
        if (count eq 1) then dataout[kk,jj] = top[igood]
        if (count eq 0) then dataout[kk,jj] = -666
    endfor
endfor

return,dataout

stop
end
