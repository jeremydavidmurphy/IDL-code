FUNCTION addnoise, spectra, noise, factor
; this function accepts an individual spectra and convolves it with a
; noise array. Both must be of the same length.

; FACTOR: This is just a multiplier to the noise, at each wavelength.

;d = readfits(spectra,h,/silent)
;n = readfits(noise,h,/silent)
d = spectra
n = noise

nn0 = n_elements(d)
s = fltarr(nn0)
for j=0,nn0-1 do begin
   noiseadd = randomn(seed,1)*(n[j])*factor
   s[j] = d[j] + noiseadd
endfor

return,s

stop
END
