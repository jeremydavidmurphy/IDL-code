; This routine is an updated version of plotsky.pro. It overplots the
; various subtracted sky values on top of the final spectra for a
; given bin. It also calculates and plots the variance...

;#####################################################################
skyfiber = 100 ;the fiber to extract from the science frame
ap = 5 ;the pixel extraction aperture
;#####################################################################

;COMPILE realignF.pro

PRO Psky_var, skylist, data

readcol,list,format='a,a,a,f,f,i',skylist,ptow,mask,wi,wdisp

test = readfits(skylist[0],/silent)
n0 = n_elements(skylist)
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*])

skyC = fltarr(n1,n0)
spectra = readfits(data,h,/silent)

for j=0,n0-1 do begin
    sky = readfits(skylist[j],h,/silent)
    ;The wavelength offset is plucked from the header
    hn = n_elements(h)
    cntr = 0
    out = 'lost'
    repeat begin
        hh = strsplit(h[cntr],' ',/extract)
        if (hh[0] eq 'WAVEZP') then out = 'found' else cntr = cntr + 1
    endrep until (out eq 'found')
    woffset = float(hh[2])
    if (woffset[j] gt 1.0) or (woffset[j] lt -1.0) then begin
        print,'Check your header! The wavelength offset is HUGE!'
        print,woffset[j]
        stop
    endif
    print,'The wavelength offset is '+woffset
    
    ;Now the files are sent into realignF...
    aligned = realignF(sky,ptow[j],mask[j],wi[j],wdisp[j],ap,woffset)

    skyC = median(aligned,dimension=2)


endfor

stop
END
