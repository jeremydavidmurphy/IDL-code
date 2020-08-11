pro bscale_fix, list

; This routine straightens out the headers coming from Vaccine which
; add a BSCALE and BZERO term. This throws off the readfits routine.

readcol,list,f='a',silent=1,files
n0 = n_elements(files)

for j=0,n0-1 do begin
    data = readfits(files[j],/silent,header)
    data = data - 2d^15
    sxdelpar,header,'BZERO'
    sxdelpar,header,'BSCALE'
    writefits,files[j],data,header
endfor

stop
END
