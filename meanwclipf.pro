; The difference here and with the meanclipF routine is that now the
; weights are included. The sigma clipping routine follows meanclipF,
; but rather than the mean being returned, the code returns a weighted
; mean of the form

;                sum (data_i * weight_i) / sum (weight_i)

; I've also stripped several pieces of the fanciness.

FUNCTION meanwclipf, data, weight, CLIPSIG=clipsig, MAXITER=maxiter, $
    CONVERGE_NUM=converge_num

compile_opt idl2
on_error,2

IF n_elements(maxiter) LT 1 THEN maxiter = 5.0
IF n_elements(clipsig) LT 1 THEN clipsig = 3.0
IF n_elements(converge_num) LT 1 THEN converge_num = 0.02

i666f = where(data eq -666,ct) ;the bad pixels are located

if (ct gt 0) then begin
    data[i666f] = !values.F_NAN 
    weight[i666f] = 0.0
endif

igoodF = where(finite(data),ct) ;index of the good values
iter = -1

REPEAT BEGIN

    goodpix = data[igoodF] ;the good data
    iter = iter + 1
    lastct = ct
    medval = median(goodpix,/even) ;JDM
    sig = stddev(goodpix,/nan)    ;JDM
    wsm = where(abs(goodpix - medval) LT clipsig * sig,ct) ;index of points kept after the clipping
    IF ct GT 0 THEN igoodF = igoodF[wsm]

ENDREP UNTIL (float(abs(ct - lastct)) / lastct LT converge_num) $
          OR (iter GT maxiter)

top = total(data[igoodF] * weight[igoodF])
bottom = total(weight[igoodF])

out = top / bottom

;print,'The number of sigma-clipping iterations was :'+strn(iter)
;print,'The weighted average is: '+strn(out)
mom = moment(data[igoodF])
mean = mom[0]
;print,'The mean is: '+strn(mean)

sigma = sqrt(mom[1])

RETURN,mean
END
