PRO readnoise

;this routine calculates the readnoise imperically, from a section of
;the CCD
readcol,'list',f='a',list
n0 = n_elements(list)

for j=0,n0-1 do begin
    file = readfits(list[j],/silent,header)
    gain = sxpar(header,'GAIN1',count=c)
    piece = file[2052:2079,100:2000]
    sd = stddev(piece)
    int = sxpar(header,'INTEGRAT',count=c2)
    rdn = sxpar(header,'RDNOISE1')
    m = mean(file)
    if (c gt 0) then begin
        print,list[j]+':  '+strn(m)+'  '+strn(gain)+'  '+strn(sd)+'  '+strn(int)
;        print,'S.D. X gain = '+strn(gain * sd)
    endif
endfor

stop
END
