PRO restack

;This works for both vp1 and vp2 and needs only the jm####.R file
;(from radplot.pro) in the calling directory. The file gets named by
;you at the end of the code...

;This code takes 1-D fits and sticks them back together into a 1D x
;247(246) fits file.
;THIS CODE REQUIRES THE ASSOCIATED jm####.R FILE TO EXIST IN THE
;CALLING DIRECTORY. It uses this information to reject the stars,dead
;fibers, etc, as determined in radplot.pro.

;It expects to find all the ##ex.fits files in the calling directory-
;As it uses the jm####.R file to determine which fibers to use, it is
;resistant to ##ex.fits files being in the directory when they
;shouldn't be.

ans=''
name=''
rv=float(0);This is the replacement value

;print,'VP1 or VP2? (Enter 1 or 2):'
;read,ans
ans = '1'
IF (ans EQ '1') then n2 = 247
;IF (ans EQ '2') then n2 = 246

jump:
print,'Enter the name of your .R file (w/o the .R):'
read,name
name = name+'.R'
file = findfile(name)
IF (file EQ '') THEN goto,jump

readcol,name,FORMAT='I,F,F',fibers,junk,cts
fnum = n_elements(fibers)

temp = fibers[bsort(fibers)]
stdfibers = intarr(n2)

name = strn(fibers[0])+'ex.fits'
test = readfits(name)
n1 = n_elements(test)

filler = intarr(n1)+rv

IF (n1 EQ 1024) THEN BEGIN
    specout=dblarr(1024,n2)
ENDIF
IF (n1 EQ 2048) THEN BEGIN
    specout=dblarr(2048,n2)
ENDIF ELSE BEGIN
    print,'The length of your ##ex.fits files is neither 1024 or 2048!?!?'
    STOP
ENDELSE

FOR j=0,n_elements(cts)-1 DO BEGIN
    IF (cts[j] EQ 0) THEN BEGIN
        print,'You have zeros in your .R file!'
        STOP
    ENDIF
ENDFOR

cntr=0
for j=0,n2-1 do begin
    if (fnum eq j-cntr) then goto,stepout
    t = temp[j-cntr]
    if (j+1 eq t) then stdfibers[j] = t
    if (j+1 ne t) then begin
        stdfibers[j] = 0
        cntr = cntr + 1
    endif
endfor
stepout:print,'Finished the loop. '+strn(n2-j)+' fibers at the top set to zero'

for j=0,n2-1 do begin
    if (stdfibers[j] ne 0) then begin
        print,'Fiber to add: '+strn(stdfibers[j])
        name = strn(stdfibers[j])+'ex.fits'
        spec = readfits(name)
        specout[*,j] = spec
    endif
    if (stdfibers[j] eq 0) then specout[*,j] = filler
endfor

;cntr=1
;FOR j=1,n2-1 DO BEGIN
;    IF (j EQ stdfibers[j-cntr]) THEN BEGIN
;        print,'Fiber to add: '+strn(stdfibers[j-cntr])
;        name = strn(stdfibers[j-cntr])+'ex.fits'
;        spec = readfits(name)
;        dout[*,j-1] = spec
;    ENDIF ELSE BEGIN
;        dout[*,j-1] = filler
;        cntr = cntr + 1
;        print,'Nothing for fiber '+strn(j) 
;        pause
;        IF (j-cntr EQ fnum) THEN BEGIN
;            FOR k = fnum,n2-1 DO BEGIN
;                dout[*,k] = filler
;                print,'fillin up the rear!'
;            ENDFOR
;        ENDIF
;    ENDELSE
;ENDFOR
print,stdfibers
print,'Name the output file (W/O THE .FITS)'
read,name
name= name+'.fits'

writefits,name,specout
STOP
END
