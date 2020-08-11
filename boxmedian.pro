FUNCTION boxmedian, data, boxsize

On_error,2

IF (n_params() EQ 0) THEN message,'Syntax -   OUT = BOXMED(data,boxsize)',/NoName

halfbox=boxsize/2
out = dblarr(n_elements(data))

FOR j=0,n_elements(data)-1 DO BEGIN
    window = data[halfbox-j:halfbox+j]
    std = window[sort(window)]
    


out = total(data)
return,out

END
