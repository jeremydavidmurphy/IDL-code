;This function accepts the actual data, the actual weight file and THE
;NAME of the mask file. It returns the weighted, collapsed data.

function wgtcollapse, data, weight, maskname, nprofile

;The quick collapse is of the form
;         collapsed = sum_over_i(pefsm_i * pefw_i)/sum_over_i(pefw_i)

;Both the data and the weight are the actual data and weight file (the
;pefsm.fits and pefw.fits files. The maskname is the name of the
;mask.dat file. The nprofile is the number of pixels in the fiber
;profile (this will always be either 5 or 7).

n1 = n_elements(data[*,0])
n2 = n_elements(data[0,*])

if (n_elements(weight[*,0]) ne n1) then begin
    print,'Your data and weight file are not the same size!'
    stop
endif
if (n_elements(weight[0,*]) ne n2) then begin
    print,'Your data and weight file are not the same size!'
    stop
endif

;To deal with the -666 flags, the corresponding weight pixels are set
;to equal zero

weight[where(data eq -666)] = 0.0

step = (nprofile-1)/2
readcol,maskname,format='i,i',y,fiber
n3 = n_elements(fiber)
n4 = nprofile
collapsed = dblarr(n1,n3)

for k=0,n3-1 do begin
    up = uint(y[k]+step)
    down = uint(y[k]-step)

; the -666 flags are maintained.    
    if (median(data[*,y[k]]) eq -666) then begin
        collapsed[*,k] = -666
        goto, jump1
    endif

    for j=0,n1-1 do begin
        t1 = total(data[j,down:up]*weight[j,down:up])
        t2 = total(weight[j,down:up])
        if (t2 eq 0) then collapsed[j,k] = -666 else begin
            collapsed[j,k] = double(t1/t2)
        endelse
    jump1:
    endfor
endfor

return,collapsed
end
