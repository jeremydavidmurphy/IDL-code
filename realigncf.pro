; This code accepts a collapsed data frame, the name of the ptow
; file, an initial wavelength, a wavelength dispersion and,
; optionally, the wavelength zeropoint offset. It completes an interpolation
; (via the INTERPOL function) and returns an aligned wavelength file.

; The inputs are the same as realignF, so both a wavelength zp
; correction and heliocentric correction can be made.
;------------------------------------------------------------------------
function realigncf, datain, ptowname, wi, wd, factor, WZP=wzp, HELIOC=helioc
;------------------------------------------------------------------------

nn1 = n_elements(datain[*,0])
nn2 = n_elements(datain[0,*])
nn1f = nn1 * factor
dataout = fltarr(nn1f,nn2) 

wo = dblarr(nn1f)
for jj=0,nn1f-1 do wo[jj] = wi + (jj*wd)/factor

wr = dblarr(nn1,nn2)
readcol,ptowname,format='d,d,d,d,d',c0,c1,c2,c3,c4

print,mean(c0)
print,mean(c1)

; the wavelength ZP correction is made
if (n_elements(wzp) ne 0) then c0 = c0 + double(wzp)

ptow = [[c0],[c1],[c2],[c3],[c4]];this is a 247x5 array, NOT a 5x247 array!

for kk=0,nn2-1 do begin
    for jj=0,nn1-1 do begin
        n = float(jj)
        wr[jj,kk] = c0[kk]+c1[kk]*n+$
          c2[kk]*n*n+c3[kk]*n*n*n+c4[kk]*n*n*n*n
    endfor
endfor

;the heliocentric correction, if any, is made
if (n_elements(helioC) ne 0) then begin
    denom = double(1.0 - (helioC / 299792.458)) ;the variable 'helioC' is now a correction factor
    wr = wr / denom
endif

for kk=0,nn2-1 do begin ;a loop over each fiber.
    temp = datain[*,kk]
    bindex = where(temp eq -666,count)
    if (count ne 0) then temp[bindex] = !values.f_nan
    
;-----------------------------------------------------------------
    dataout[*,kk] = interpol(temp,wr[*,kk],wo)
;-----------------------------------------------------------------
endfor

;the nan's and inf's are changed back to -666 flags.    
p = finite(dataout)
ibad = where(p eq 0,count)
if (count gt 0) then dataout[ibad] = -666.0

return,dataout

END
