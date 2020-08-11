; this routine accepts a 3D array and returns the median (along the
; 3rd dimension) as the output. It is good for median combining a
; series of frames of the same image.

function med3d, array

n1 = n_elements(array[*,0,0])
n2 = n_elements(array[0,*,0])
n3 = n_elements(array[0,0,*])

out = dblarr(n1,n2)
for j=0,n1-1 do begin
    for k=0,n2-1 do begin
        out[j,k] = median(array[j,k,*],/even)
    endfor
endfor

return,out

end
