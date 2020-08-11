; This function finds the minimum X2 values for a given parameter. It
; returns the minimum X2 value at each of the given parameters.

FUNCTION findx2floorf, param, x2

; PARAM: The parameter values
; X2: The associated chi-squared values

n0 = n_elements(param)
iparam = bsort(param)
p = param[iparam] ;parameters are sorted from low to high

x2s = x2[iparam] ;the chi-squared values are also sorted, based on the low-to-high parameter values

paramX = [0.0]
x2Y = [0.0]
i = 0
repeat begin
    p1 = p[i]
    p2 = p[i+1]
    if (p2 gt p1) then begin
        ii = where(param eq p1)
        paramX = [mlx0,p1]
        x2Y = [x2Y,min(x2[ii])]
    endif
    i = i + 1
endrep until (p2 ge max(param))
ii = where(param eq p2)
paramX = [paramX,p2]
x2Y = [x2Y,min(param[ii])]
paramX = paramX[1:*]
x2Y = x2Y[1:*]

out = [[paramX],[x2Y]]
return,out

stop
END
