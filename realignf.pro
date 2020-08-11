;THE OLD VERSION THAT RAN ON COLLAPSED DATA IS NOW SAVED AS
;REALIGNC.PRO. THE NEW VERSION ALLOWS FOR UNCOLLAPSED DATA. IT ALSO
;EXPECTS THE dotR FILES TO ALREADY HAVE BEEN APPLIED (AS IT IS IN
;PIPE2.PRO) SO THAT THE REJECTED FIBERS HAVE BEEN MASKED WITH THE -666
;FLAG.

;The returned file is the same size as the input file.
;----------------------------------------------------------------------------
function realignf, datain, ptowframename, maskframename, wi, wdisp, $
                   aperture, Factor=factor, WZP=wzp, helioC=helioC, VAC=vac
;compile_opt idl2
;on_error,2
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
;GETS ADDED to the zeroth order term in the wavelength solution. (I
;double-checked this on July 8, 2011- it's the correct direction.)
;HELIOC: The heliocentric correction. The number that Vaccine returns
;is THE EARTH'S VELOCITY IN THE DIRECTION OF THE OBJECT YOU'RE
;OBSERVING. So, this means that a positive vhelio value has the earth
;moving towards an object = a blue shift. To remove this you use the
;following equation

;wavelength_rest = wavelength_observed / ( 1 - vhelio/c)
;(this equation was also double-checked on July 11, 2011 and found to
;be correct)

;The equation for the wavelength realignment is
;
;             waveout(x) = WI + WDISP * x
;
;where 'x' is the pixel position.

; TO TURN ON WZP AND/OR HELIOC THEN SET THEM EQUAL TO A NUMBER OTHER
; THAN ZERO (NOT A STRING).

;FACTOR: Added on Dec 28th, 2010 to allow for an interpolation
;factor. This allows for a super-sampling of the data by a given
;factor.

;THE NEW SOLUTION IS DETERMINED BY THE IDL ROUTINE INTERPOL.PRO. 

;AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4)
;The above equation is the one SDSS uses and is the IAU standard as
;laid out in Morton (1991, ApJS, 77, 119)
;This equation can be simplified to n = 1.000279 with air wavelengths
;being shorter (bluer) than vacuum wavelengths. Another source to look
;is the SDSS website at
;http://www.sdss.org/dr7/products/spectra/vacwavelength.html
;Josh uses n = 1.00022, based on a Filipenko paper that did the
;calculation for different latitudes and conditions. This is his
;estimate of the average conditions at McDonald.

;VAC: set to 'yes' if you want to convert to vacuum wavelengths.

step = floor((aperture-1)/2.0)

readcol,silent=1,maskframename,format='i,x',yindex
yindex = yindex-1

nn1 = n_elements(datain[*,0])
nn2 = n_elements(datain[0,*])

if (n_elements(factor) gt 0) then begin
    nn1f = nn1 * factor 
    print,'Increasing original array size of '+strn(nn1)+' to '+strn(nn1f)
endif else nn1f = nn1

;The bad fibers are located. This assumes that fewer than 1/2 of a
;given row of data is -666 flags as it uses the median of a row to
;determine which fibers have been masked.
;tm = median(datain,dimension=1)
;badfibers = where(tm eq -666)

wr = dblarr(nn1,nn2)
wo = dblarr(nn1f)
dataout = dblarr(nn1f,nn2)

readcol,silent=1,ptowframename,format='d,d,d,d,d',c0,c1,c2,c3,c4

; The zeropoint correction, if any, is made
if (n_elements(wzp) ne 0) then c0 = c0 + double(wzp)

ptow = [[c0],[c1],[c2],[c3],[c4]];this is a #_of_fibers X 5 array

;the wavelength (real) array is created from the ptow.dat files
for k=0,n_elements(c0)-1 do begin
    for j=0,nn1-1 do begin
        n = float(j)
        wr[j,yindex[k]-step:yindex[k]+step] = $
          c0[k] + c1[k]*n + c2[k]*n*n + c3[k]*n*n*n + c4[k]*n*n*n*n
    endfor
endfor

;the heliocentric correction, if any, is made
if (n_elements(helioC) ne 0) then begin
    denom = double(1.0 - (helioC / 299792.458))
    wr = wr / denom
endif

if (n_elements(VAC) ne 0) then begin
    if (VAC eq 'yes') then begin
        wr = wr * 1.000279
        print,'Converting the wavelength to VACUUM...'
    endif
endif

;the wavelength out array is created
for j=0,nn1f-1 do wo[j] = wi + (j*wdisp)/factor

;all rows that are all -666 flags (bad fibers) or zeros (gaps) are
;returned as the original value.
for k=0,nn2-1 do begin ;a loop through each row of data
    test = median(datain[*,k])
    if (test eq -666.0) then begin
        dataout[*,k] = -666.0 ;fibers that were rejected during pipe1.pro 
        goto, jump1
    endif
    if (test eq 0) then begin
        dataout[*,k] = 0 ;dead fibers in the IFU
        goto, jump1
    endif

;any remaining -666 flags from cosmics are set to 'NAN' BEFORE going
;through the interpolation. Interpol.pro handles these as missing data
;and returns 'nan'. This gets around the issue of interpolating
;between a real value and a -666 flag, which would return some bad and
;likely (BUT NOT ALWAYS) negative value. The 'nan' flags are then set
;back to -666 flags and are returned to the pipeline.

    temp = datain[*,k]
    bindex = where(temp eq -666.0,count)
    if (count gt 0) then temp[bindex] = !Values.F_NAN
;    print,'Interpolating row #'+strn(k+1)

;-----------------------------------------------------------------
    dataout[*,k] = interpol(temp,wr[*,k],wo)
;-----------------------------------------------------------------

    jump1:
endfor

p = finite(dataout)
ibad = where(p eq 0,count)
if (count gt 0) then dataout[ibad] = -666.0

return,dataout
end
