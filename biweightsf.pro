FUNCTION biweightSF, array
ON_ERROR,2

; The "S" stands for single. This code takes in a set of numbers and
; returns the biweight and biweight uncertainty.

; The modification to this code (as compared to biweightF.pro) is just
; to allow a 1-D array of 
; numbers to be read in (as opposed to spectra) and biweighted. If an
; array larger than a 1-D 
; array is sent in, the array will be made into a 1-D array and the
; biweight returned

; MODIFIED ON MAY 12, 2011: This code now returns the uncertainty of
; the value, as per the 1990 Beers paper.

;THE CONCENTRATION PARAMETER****************************************
;The HIGHER the number the WEAKER the rejection
c =  6.0 ;The concentration (tension) parameter.
cE = 9.0 ;The uncertainty concentration parameteprint,'Running biweightSF...'
speed = 10e-6 ;10e-4 is very fast, 10e-6 is slow, and slightly more 
              ;robust for arrays with small numbers
;*******************************************************************

s = size(array)
if (s[0] gt 1) then array = array[bsort(array)]

biwt = dblarr(2)

;any -666 flags are converted to NaN's and thus rejected by the median
ibadF = where(array eq -666,countF)
if (countF gt 0) then array[ibadF] = !values.F_NaN

temp = array
pp = finite(temp)
ii = where(pp eq 1,countF)

if(countF gt 1) then begin
    temp = temp[ii]
    nn4 = n_elements(temp)
    nn4E = float(nn4)
    M = median(temp,/even,/double)
    goto,jump1
endif
if (countF eq 1) then begin
    biwt[0] = temp[ii]
    biwt[1] = 0.0
    goto,jumpend
endif
if (countF eq 0) then begin
    biwt[0] = temp[0]
    biwt[1] = temp[0]
    print,'NaNs for every element in the array!!!'
    goto,jumpend
endif
    
jump1:
;------------------------------------------------------------------
biarr1 = dblarr(nn4)
biarr2 = dblarr(nn4)
biarr1E = dblarr(nn4)
biarr2E = dblarr(nn4)
madarr = dblarr(nn4)
ui = dblarr(nn4)
uiE = dblarr(nn4)
cntrF = 0
delta = 1.0
;--------------------------------------------------------------------------------------
;THE BIWEIGHT
;--------------------------------------------------------------------------------------
while (abs(delta) gt speed and cntrF lt 20) do begin
    madarr = abs(temp - M) ;the MAD array, BEFORE the median is taken
    MAD = median(madarr,/even,/double)
    if (MAD eq 0.0) then MAD = 0.0000000001

    kk = 0.0
    repeat begin
        ui[kk] = (temp[kk]-M)/(c*MAD)
        t = abs(ui[kk])
        if (t le 1.0) then begin ;if the point receives a weight
            biarr2[kk] = (1-ui[kk]^2)^2
            biarr1[kk] = (temp[kk]-M)*biarr2[kk]
            biarr1E[kk] = ((temp[kk]-M)^2) * ((1-(uiE[kk])^2)^4)
            biarr2E[kk] = (1-(uiE[kk]^2)) * (1-(5*(uiE[kk]^2)))
        endif else begin        ;if the point is rejected
            biarr2[kk] = 0.0
            biarr1[kk] = 0.0
            biarr1E[kk] = 0.0
            biarr2E[kk] = 0.0
        endelse
        kk = kk + 1.0
    endrep until (kk eq nn4)
    
    t1 = total(biarr1)
    t2 = total(biarr2)
    t1E = sqrt(total(biarr1E)) * sqrt(nn4E)
    t2E = abs(total(biarr2E))

    if (t2 eq 0.0) then begin
        biwt[0] = -666.0
        biwt[1] = -666.0
        goto,jumpend
    endif else delta = t1/t2
    M = M + delta
    cntrF = cntrF + 1
    
endwhile

biwt[0] = M
biwt[1] = double(t1E/t2E)

print,'The number of biweight iterations was '+strn(cntrF)
print,'Delta equals '+strn(delta)

jumpend:

RETURN,biwt

END
