; This function accepts the center of the galaxy, and the fiber
; RA and Dec, then calcuates the RA and Dec offsets for each
; fiber. The point of this routine is that it accepts the coordinate
; in the standard RA and Dec format (e.g. 23:12:34.8 -14:09:34) AS
; SCRIPTS, then does the conversion.

; RAoff and Decoff are offsets, in arcseconds, that you can add to
; adjust the output values. THEY ARE SUBTRACTED FROM THE FINAL dRA and
; dDec values that are returned.

; The returned values are (in arcseconds):
; galRA   galDec   fibRA   fibDec   dRA   dDec   radius

FUNCTION offsetF, galRA, galDec, fiberRA, fiberDec, $
                 RAoff = raoff, Decoff = decoff

RA = strn(galRA) & Dec = strn(galDec)
RA = strsplit(RA,':',/EXTRACT)
Dec = strsplit(Dec,':',/EXTRACT)

out1 = double((RA[0]+RA[1]/60.+RA[2]/3600.)*15)
;The + or - sign of the declination is accounted for
if (Dec[0] lt 0) then out2 = double(Dec[0]-Dec[1]/60.-Dec[2]/3600.)
if (Dec[0] ge 0) then out2 = double(Dec[0]+Dec[1]/60.+Dec[2]/3600.)

fRA = strsplit(fiberRA,':',/EXTRACT)
fDec = strsplit(fiberDec,':',/EXTRACT)

out3 = double((fRA[0]+fRA[1]/60.+fRA[2]/3600.)*15.)

if (fDec[0] lt 0) then out4 = double(fDec[0]-fDec[1]/60.-fDec[2]/3600.)
if (fDec[0] gt 0) then out4 = double(fDec[0]+fDec[1]/60.+fDec[2]/3600.)


out5 = double((out3 - out1) * 3600.0)
out6 = double((out4 - out2) * 3600.0)

if (n_elements(raoff) ne 0) then begin
    raoff = double(raoff)
    out5 = out5 - raoff
endif
if (n_elements(decoff) ne 0) then begin
    decoff = double(decoff)
    out6 = out6 - decoff
endif

out7 = double(sqrt(out5^2 + out6^2))

out = [out1,out2,out3,out4,out5,out6,out7]
return,out


END
