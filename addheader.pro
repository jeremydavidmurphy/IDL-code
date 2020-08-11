PRO addheader, list

; This routine adds header elements in fits files. 

; LIST: the list of fits files you want headers added to...

readcol,list,f='a',silent=1,files
print,files
n0 = n_elements(files)

for j=0,n0-1 do begin
    data = readfits(files[j],header,/silent)
;    sxaddpar,header,'COMMENT','ROW1: Wavelength'
;    sxaddpar,header,'COMMENT','ROW2: Signal-2-Noise'
;    sxaddpar,header,'COMMENT','ROW3: Instrumental Sigma (km/sec)'
;    sxaddpar,header,'COMMENT','ROW4: Instrumental Resolution (FWHM in A)'
;    sxaddpar,header,'COMMENT','ROW5: Normalized flux'
    sxaddpar,header,'RDNOISE1',3.7,' a fudge to get Vaccine to work. (was 4.2)'
    sxaddpar,header,'INTEGRAT',2,' a fudge to get Vaccine to work. (was 3)'
    writefits,files[j],data,header
endfor

stop
END
