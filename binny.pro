; OBSOLETE

PRO binny

;This routine gets applied AFTER the frames have been combined in
;IRAF, leaving a final fits file that's ????x247
;
;The routine lumps the data into it's bins, then completes a weighted
;combination.

;This code expects the dead fibers to be present, with a value of
;0. They are included in the combination, yet combined with a weight
;of 0.

name=''
ans=''
print,'Enter the name of the fits file to bin:'
read,name
name = name+'.fits'


datain = readfits(name)
denom = dblarr(n_elements(datain[*,0]))
files = findfile('bin????',count=num)
spectra = dblarr(n_elements(datain[*,0]),num)

FOR j=1,num DO BEGIN
    name = files[j-1]
    readcol,name,Format='I',list
    n = n_elements(list)
    chunk = dblarr(n_elements(datain[*,0]),n)
    weight = chunk
    wtdata = chunk
    FOR k=0,n-1 DO BEGIN ;This puts all the fibers for the given bin into 1 array
        chunk[*,k] = datain[*,list[k]-1]
    ENDFOR

    FOR k=0,2047 DO denom[k] = total(chunk[k,*])

    FOR k=0,n-1 DO BEGIN
        weight[*,k] = chunk[*,k]/denom
        wtdata[*,k] = weight[*,k]*chunk[*,k]
    ENDFOR
    FOR k=0,2047 DO spectra[k,j-1] = total(wtdata[k,*])
ENDFOR

print,'Name the output file (W/O THE .FITS):'
read,name
name = name+'.fits'

writefits,name,spectra

STOP
END
