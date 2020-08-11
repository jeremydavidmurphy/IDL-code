PRO sumbins, dataname

;This code reads in a data frame (dataname) and a list of bin
;names. It then completes a STRAIGHT SUM. It uses the names in the bin list
;(expected to be of the form ####.bin and writes the new spectra to
;####.fits)

binlist=''
ans=''
print,'Enter the name of your bin list:'
read,binlist

datain = readfits(dataname)
n1 = n_elements(datain[*,0]);wavelength


readcol,binlist,FORMAT='A',fbin
n2 = n_elements(fbin)
databin = dblarr(n1,n2) ;
outnames = strarr(n2)
for j=0,n2-1 do begin
    nm = strsplit(fbin[j],'.',/extract)
    outnames[j] = nm[0]
endfor

FOR j=0,n2-1 DO BEGIN ;the 'bin' loop is run

    readcol,fbin[j],FORMAT='I',list
    list = list-1
    n3 = n_elements(list)

    IF (n3 EQ 1) THEN BEGIN;a bin size of ONE
        ind = list
        databin[*,j] = datain[*,ind]
    ENDIF
    IF (n3 GT 1) THEN BEGIN
        chunk = dblarr(n1,n3)
;        wgts = dblarr(n3)
        FOR k=0,n3-1 DO BEGIN
            ind = list[k]
            chunk[*,k] = datain[*,ind]
;            wgts[k] = median
        ENDFOR
        databin[*,j] = total(chunk,2)
    ENDIF

ENDFOR

ans=''
print,'Lump them together or output as individual fits files?'
print,'(type "lump" or "one")'
read,ans

if (ans eq 'lump') then begin
    name = ''
    print,'Name the output fits file:'
    read,name
    writefits,name,databin
endif

if (ans eq 'one') then begin
    for k=0,n2-1 do begin
        writefits,outnames[k]+'.fits',databin[*,k]
        print,'File written as '+outnames[k]+'.fits...'
    endfor
endif
    
STOP
END
