PRO restack_vp1

;This code takes 1-D fits and sticks them back together into a 1D x
;247 fits file.
;THIS CODE REQUIRES THE ASSOCIATED jm####.R FILE TO EXIST IN THE
;CALLING DIRECTORY. It uses this information to reject the stars, etc,
;as determined in radplot.pro.

;It expects to find all the ##ex.fits files in the calling directory-
;As it uses the jm####.R file to determine which fibers to use, it is
;resistant to ##ex.fits files being in the directory when they
;shouldn't be.

ans=''
name=''
rv=float(0);This is the replacement value
n2=0

files = findfile('*ex.fits',count=num)

dead = [16, 141, 200, 215, 233, 234]-1

test = readfits(files[0])
n1 = n_elements(test)

;Just a check to confirm your data is the right size
IF (n1 EQ 1024) THEN BEGIN
    dout=dblarr(1024,247)
ENDIF
IF (n1 EQ 2048) THEN BEGIN
    dout=dblarr(2048,247)
ENDIF ELSE BEGIN
    print,'The length of your ##ex.fits files is neither 1024 or 2048!?!?'
    STOP
ENDELSE

print,'Enter the name of your .R file (w/o the .R):'
read,name
name = name+'.R'

readcol,name,FORMAT='I,F,F',fibers,junk,cts
fnum = n_elements(fibers)

FOR j=0,n_elements(cts)-1 DO BEGIN
    IF (cts[j] EQ 0) THEN print,'You have zeros in your .R file!'
    
ENDFOR

stdfibers = fibers[bsort(fibers)]
print,stdfibers
pause

cntr=0
filler = intarr(n1)+rv
FOR j=0,fnum-1 DO BEGIN
    IF (j EQ dead[cntr]) THEN BEGIN
        dout[*,j] = filler
        cntr = cntr+1
        IF (cntr EQ n_elements(dead)) THEN cntr=0
    ENDIF ELSE BEGIN
        name = strn(j+1)+f
        spec = readfits(name)
        dout[*,j] = spec
    ENDELSE
ENDFOR

;THIS PART IS STILL BROKEN... 

;cntr=0
;dead=dead+1
;IF (ans EQ 'n') THEN BEGIN
;    FOR j=1,247 DO BEGIN
;        IF (j NE dead[cntr]) THEN BEGIN
;            name = strn(j)+f
;            print,name
;            spec = readfits(name)
;            dout[*,j-1-cntr] = spec
;        ENDIF ELSE BEGIN
;            cntr = cntr+1
;            IF (cntr EQ n_elements(dead)) THEN cntr=0
;        ENDELSE
;    ENDFOR
;ENDIF

;This next piece reads in the ###ex.fits file names and sorts them
;into numerical order

;ind = intarr(241)
;FOR j=0,240 DO BEGIN
;    temp = strsplit(files[j],suf,/extract)
;    ind[j] = uint(temp[0])
;    ind[j] = ind[j]-1
;ENDFOR
;ind = ind[sort(ind)]
;FOR j=0,240 DO files[j] = files[ind[j]]; The files list gets sorted into order

;dout = dblarr(2048,247)

;print,'Reinsert the dead fibers?'
;read,ans
;IF (ans EQ 'y') THEN BEGIN
;cntr=0
;print,'At what value?'
;read,rv
;filler = intarr(2048)+rv
;    FOR j=1,247 DO BEGIN
;        temp = strsplit(files[j-1-cntr],suf,/extract)
;        temp = uint(temp[0])
;        print,temp & pause
;        IF (temp EQ j) THEN BEGIN
;            dout[*,j-1] = readfits(files[j-1-cntr])
;        ENDIF ELSE BEGIN
;            dout[*,j-1] = filler
;            cntr = cntr + 1
;        ENDELSE
;    ENDFOR
;ENDIF

print,'Are there fibers to reject?'
read,ans
IF (ans EQ 'y') THEN BEGIN
    print,'Enter the number of fibers to reject:'
    read,n2
    dump=intarr(n2)
    print,'Enter the actual fibers to reject (EX: 23 46 47 129):'
    read,dump
    dump = dump-1
    FOR j=0,n2-1 DO BEGIN
        cntr = dump[j]
        dout[*,cntr] = filler 
    ENDFOR
ENDIF ELSE print,'No fibers will be rejected.'



;These two methods are similar.

;FOR j=0,num-1 DO BEGIN
;    IF (j NE dump[cnt]) THEN BEGIN
;        name = strn(j)+f
;        data[*,j] = readfits(name)
;    ENDIF ELSE BEGIN
;        data[*,j] = dead
;        cnt = cnt+1
;    ENDELSE
;ENDFOR

;FOR j=1,num DO BEGIN
;    name = strn(j)+f
;    data[*,j-1] = readfits(name)
;ENDFOR

;IF (n_elements(dump) NE 0) THEN BEGIN
;    FOR j=0,n_elements(dump)-1 DO BEGIN
;        data[*,dump[j]-1] = dead
;    ENDFOR
;ENDIF


print,'Name the output file (W/O THE .FITS)'
read,name
name= name+'.fits'

writefits,name,dout
STOP
END
