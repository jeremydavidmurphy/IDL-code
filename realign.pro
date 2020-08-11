;THE OLD VERSION THAT RAN ON COLLAPSED DATA IS NOW SAVED AS
;REALIGNC.PRO. THE NEW VERSION ALLOWS FOR UNCOLLAPSED DATA. IT ALSO
;EXPECTS THE dotR FILES TO ALREADY HAVE BEEN APPLIED (AS IT IS IN
;PIPE2.PRO) SO THAT THE REJECTED FIBERS HAVE BEEN MASKED WITH THE -666 FLAG.
;--------------------------------------------------------------------------
pro realign, dataname, ptowframename, maskframename, wi, wdisp, aperture, $
             WZP=wzp, OUTPUT=KorB
;--------------------------------------------------------------------------

; MODIFIED (APRIL 8, 2009): This code now accepts another entry. This
; comes from the header and makes a correction to the ptow.dat values
; as a linear shift based on a sky line as determined by Vaccine and
; written into the header.

; DATANAME: The name of the data, W/O the .fits. This will be a MxN array.
; PTOWFRAMENAME: The name of the ptow.dat frame. REALIGN.PRO then looks
; for this file in the calling directory
; MASKFRAMENAME: Same as above.
; WI: The initial wavelength for the wavelength linearization as given
; in the equation below
; WDISP: The wavelength dispersion.
; APERTURE: The size of the extraction window. SET TO ZERO FOR
; COLLAPSED DATA.

; KEYWORD CALLS:
; WZP: This keyword can be used in a couple ways. If WZP='find' then the
; routine will extract the wavelength zero point from the header of
; the file and ADD it to the zeroth term of the ptow.dat files.If WZP
; = ### then that number will be used as the zero point offset.

; OUT: This determines whether the spectra are returned as one file
; (OUT='keep'), with the original name being subscripted with an
; "A". Or, if OUT='break' then each fiber is written out individually,
; as dataname_###.fits. IF THE OUT CALLWORD IS NOT USED, THEN JUST THE
; ONE FILE IS WRITTEN OUT AND SUBSCRIPTED WITH AN "A" JUST AS WITH OUT
; = "keep"

; The equation for the wavelength realignment is
;
;             waveout(x) = WI + WDISP*x
;
; where 'x' is the pixel position.

; THE NEW SOLUTION IS DETERMINED BY THE IDL ROUTINE INTERPOL.PRO. 

;MODIFIED Nov 15, 2009: The "break" statement is now working
;properly. The code now finds the corresponding weight file and uses
;this to make a weighted collapse of the data, over the size of the
;aperture. AS NEGATIVE VALUES ARE UNPHYSICAL, THEY ARE SET TO ZERO
;BEFORE THE COLLAPSE.

if (aperture ne 0) then begin
    step = floor((aperture-1)/2.)
    col = 'col'
    readcol,silent=1,maskframename,format='i,x',yindex
    yindex = yindex-1
endif else begin
    step = 0
    col = 'uncol'
endelse

data = readfits(dataname+'.fits',/silent,header)
if (n_elements(wzp) ne 0) then begin
    if (wzp eq 'find') then begin
        wavezp = 'on'
        cntr3 = 0
        repeat begin
            temp = strsplit(header[cntr3],' ',/extract)
            ans = temp[0]
            cntr3 = cntr3 + 1
        endrep until (ans eq 'WAVEZP')
        wzp = double(temp[2])
        print,'The wavelength zeropoint is: '+strn(wzp)
        print,'If this is not a number, correct your code...'
        pause
    endif
endif

nn1 = n_elements(data[*,0])
nn2 = n_elements(data[0,*])
nn3 = n_elements(yindex)

;The bad fibers are located. This assumes that fewer than 1/2 of a
;given row of data is -666 flags as it uses the median of a row to
;determine which fibers have been masked.
tm = median(data,dimension=1)
badfibers = where(tm eq -666)

wr = dblarr(nn1,nn2) ;the "real" wavelength file
wo = dblarr(nn1) ;the "output" wavelength file
dataout = dblarr(nn1,nn2)

readcol,silent=1,ptowframename,format='d,d,d,d,d',c0,c1,c2,c3,c4

if (n_elements(wzp) ne 0) then c0 = c0+wzp

ptow = [[c0],[c1],[c2],[c3],[c4]];this is a 247x5 array, NOT a 5x247 array!

;the wavelength (real) array is created from the ptow.dat files
for k=0,n_elements(c0)-1 do begin
    if (col eq 'col') then begin
        for j=0,nn1-1 do begin
            wr[j,yindex[k]-step:yindex[k]+step] = $
              c0[k]+c1[k]*j+c2[k]*j*j+c3[k]*j*j*j+c4[k]*j*j*j*j
        endfor
    endif
    if (col eq 'uncol') then begin
        for j=0,nn1-1 do begin
            wr[j,k] = c0[k]+c1[k]*j+$
              c2[k]*j*j+c3[k]*j*j*j+c4[k]*j*j*j*j
        endfor
    endif
endfor

;the wavelength out array is created
for j=0,nn1-1 do wo[j] = wi + (j*wdisp)

;all rows that are all -666 flags (bad fibers) or zeros (gaps) are
;returned as the original value.
for k=0,nn2-1 do begin
    test = median(data[*,k])
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

    temp = data[*,k]
    bindex = where(temp eq -666,count)
    if (count ne 0) then temp[bindex] = !values.f_nan
    dataout[*,k] = interpol(temp,wr[*,k],wo)
;the nan's and inf's are changed back to -666 flags.    
    for j=0,nn1-1 do begin
        p = finite(dataout[j,k])
        if (p eq 0) then dataout[j,k] = -666
    endfor

    jump1:
endfor

if (n_elements(KorB) eq 0) then begin
    writefits,dataname+'_A.fits',dataout
    goto,jump2
endif

if (KorB eq 'keep') then begin
    writefits,dataname+'_A.fits',dataout
    goto,jump2
endif

;you need to use the mask file here and output collapsed files here...
if (KorB eq 'break') then begin
    temp = strsplit(dataname,'_',/extract)
    weights = readfits(temp[0]+'_pefw.fits')
    for j=0,nn3-1 do begin
        piece = dataout[*,yindex[j]-step:yindex[j]+step]
        piece[where(piece lt 0.0)] = 0.0
        piecew =  weights[*,yindex[j]-step:yindex[j]+step]
        if (median(piece) ne 0 or median(piece) ne -666) then begin
            temp = strsplit(dataname,'_',/extract)
            name = temp[0]+'_'+strn(j+1)+'.fits'
            pieceout = piece * piecew
            pieceC = total(pieceout,2)
            writefits,name,pieceC
            print,'fiber # '+strn(j+1)+' was written out...'
        endif
    endfor
endif
jump2:
stop
end
