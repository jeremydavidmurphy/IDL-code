; MODIFIED FROM REALIGNF.PRO ON DEC 26, 2010. This routine is altered
; from the original realignF.pro in that it super-samples the
; interpolation by a factor of 3 (it's actually a parameter that's
; hardcoded just below. So, insert a 2000 element array and
; it will return a 6000 element array.

;The returned file is the same size as the input file.
;----------------------------------------------------------------------------
function realignF, datain, ptowframename, maskframename, wi, wdisp, $
                   aperture, WZP=wzp, helioC=helioC
compile_opt idl2
on_error,2

factor = 3.0 ;change this to change the super-sampling factor.
;----------------------------------------------------------------------------

; MODIFIED (APRIL 8, 2009): This code now offers another entry. This
; comes from the header and makes a correction to the ptow.dat values
; as a linear shift based on a sky line as determined by Vaccine and
; written into the header.

; MODIFIED (AUG 3, 2010): The heliocentric correction is now an
; option. If used, this value will get ADDED to the first term in the
; wavelength solution.

;DATAIN: The actual data to be realigned. This will be a MxN array.
;PTOWFRAMENAME: The name of the ptow.dat frame. REALIGN.PRO then looks
;for this file in the calling directory
;MASKFRAMENAME: Same as above.
;WI: The initial wavelength for the wavelength linearization as given
;in the equation below
;WDISP: The wavelength dispersion.
;APERTURE: The number of pixels in the Y-direction per fiber.
;WZP: The number taken from the fits header labeled WZP. This number
;GETS ADDED to the zeroth order term in the wavelength solution.
;HELIOC: The heliocentric correction. This is in km/sec and so, if
;used, get's ADDED to each fiber's wavelength solution BEFORE
;interpolation. It get's handled via

; helio_correction(wavelength) = (helioC / c) * wavelength

;The equation for the wavelength realignment is
;
;             waveout(x) = WI + WDISP*x
;
;where 'x' is the pixel position.

;THE NEW SOLUTION IS DETERMINED BY THE IDL ROUTINE INTERPOL.PRO. 

step = floor((aperture-1)/2.0)

readcol,silent=1,maskframename,format='i,x',yindex
yindex = yindex-1

nn1 = n_elements(datain[*,0])
nn2 = n_elements(datain[0,*])

;The bad fibers are located. This assumes that fewer than 1/2 of a
;given row of data is -666 flags as it uses the median of a row to
;determine which fibers have been masked.
;tm = median(datain,dimension=1)
;badfibers = where(tm eq -666)

wr = dblarr(nn1,nn2)
wo = dblarr(nn1)
dataout = dblarr(nn1,nn2)

readcol,silent=1,ptowframename,format='d,d,d,d,d',c0,c1,c2,c3,c4

; The zeropoint correction, if any, is made
if (n_elements(wzp) ne 0) then begin
    c0 = c0 + wzp
endif

ptow = [[c0],[c1],[c2],[c3],[c4]];this is a 247x5 array, NOT a 5x247 array!

;the wavelength (real) array is created from the ptow.dat files
for k=0,n_elements(c0)-1 do begin
    for j=0,nn1-1 do begin
        n = float(j)
;        wr[j,yindex[k]-step:yindex[k]+step] = $
;          c0[k] + c1[k]*n + c2[k]*n^2 + c3[k]*n^3 + c4[k]*n^4
        wr[j,yindex[k]-step:yindex[k]+step] = $
          c0[k] + c1[k]*n + c2[k]*n*n + c3[k]*n*n*n + c4[k]*n*n*n*n
    endfor
endfor

;the heliocentric correction, if any, is made
if (n_elements(helioC) ne 0) then begin
    helioC = 1.0 - (helioC / 299792.458) ;the variable 'helioC' is now a correction factor
    wr = wr / helioC
endif

;the wavelength out array is created
;for j=0,nn1-1 do wo[j] = wi + (j*wdisp)
wo = dindgen(nn1*factor) * (wdisp/factor) + wi

;all rows that are all -666 flags (bad fibers) or zeros (gaps) are
;returned as the original value.
for k=0,nn2-1 do begin
    test = median(datain[*,k])
    if (test eq -666) then begin
        dataout[*,k] = -666
        goto, jump1
    endif
    if (test eq 0) then begin
        dataout[*,k] = 0
        goto, jump1
    endif

;any remaining -666 flags from cosmics are set to 'inf' BEFORE going
;through the interpolation. Interpol.pro handles these as missing data
;and returns 'nan'. This gets around the issue of interpolating
;between a real value and a -666 flag, which would return some bad and
;likely (BUT NOT ALWAYS) negative value. The 'nan' flags are then set
;back to -666 flags and are returned to the pipeline.

    temp = datain[*,k]
    bindex = where(temp eq -666,count)
    if (count gt 0) then temp[bindex] = !Values.F_NAN

;-----------------------------------------------------------------
    dataout[*,k] = interpol(temp,wr[*,k],wo)
;-----------------------------------------------------------------

p = finite(dataout)
ibad = where(p eq 0,count)
if (count gt 0) then dataout[ibad] = -666.0

return,dataout
end
