FUNCTION medbox, spec, boxsize, medorave
ON_ERROR,2

;This is a median boxcar routine, where the median at the start and
;end are taken from a increasing[decreasing] box until the boxsize is reached

;MY FIRST FUNCTION
;Accepts a 1-D spectra and returns a median boxcar smoothed spectra
;If no boxsize is input it sets a default boxsize of 21

;IF (n_elements(boxsize) EQ 0) THEN boxsize=21
;The switch between median or average boxcar
IF (medorave EQ 1) THEN select=1
IF (medorave EQ 2) THEN select=2

halfbox=boxsize/2
out = dblarr(n_elements(spec))
n = n_elements(out)

;Returns a median boxcar
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF (select EQ 1) THEN BEGIN
    FOR j=0,halfbox-1 DO out[j] = median(spec[0:j],/EVEN)
    FOR j=halfbox,n-halfbox-1 DO out[j] = median(spec[j-halfbox:j+halfbox],/EVEN)
    FOR j=n-halfbox,n-1 DO out[j] = median(spec[j:n-1],/EVEN)
ENDIF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Returns an average boxcar
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF (select EQ 2) THEN BEGIN
    FOR j=0,halfbox-1 DO BEGIN
        nn = n_elements(spec[0:j])
        out[j] = total(spec[0:j])/nn
    ENDFOR
    FOR j=halfbox,n-halfbox-1 DO BEGIN
        nn = n_elements(spec[j-halfbox:j+halfbox])
        out[j] = total(spec[j-halfbox:j+halfbox])/nn
    ENDFOR
    FOR j=n-halfbox,n-1 DO BEGIN
        nn = n_elements(spec[j:n-1])
        out[j] = total(spec[j:n-1])/nn
    ENDFOR
ENDIF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

RETURN,out
END
