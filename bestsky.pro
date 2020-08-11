PRO bestsky, list, flip

; This routine is a 2-in-1.
; If "flip" is set to 'one' then the first routine is run. This uses
; the pefy, pe, and flat files to generate a bunch of temporary,
; collapsed pefy files, based on different multiples of the level of
; sky-subtraction. THIS ROUTINE HAS NOT BEEN WRITTEN YET.

; If "flip" is set to 'two' then this routine reads in a number of
; already reduced pefsm.fits files

; This routine is used to do a weighted biweight sum over an entire 
; frame. The purpose of this is to then fit with a fixed template star,
; and assume the minimum RMS of that fit is the best level of 
; sky-subtraction.

; PEFS = (PE/FLAT) - PEFY

; The idea here is to estimate the best sky-subtraction BEFORE using
; Vaccine to subtract the sky. The way this will be done in future runs
; is to use Vaccine to generate the pefy.fits files for AN INDIVIDUAL
; SKY FRAME. Then, apply the above equation to estimate the best sky
; with multiples of the two bracketing sky frames. Once the range is
; pinned down in this manner, Vaccine is used to generate the final
; data frames.

; The lists are of the same format as those fed into pipe2.pro.
; EX:
; jm0300cc_   b   _jan08_n5.dat   jm0300.R



if (flip eq 'one') then begin

endif

if (flip eq 'two') then begin
    
;****************************************************************************
    wi = 3530.0
    wdisp = 1.12
    step = 3
; STEP = (NPROFILE-1)/2  set to 3 for 7 pixel and to 2 for 5 pixel extraction
;****************************************************************************
    readcol,silent=1,list,format='a,x,a,a',files,dat,dotR,silent=1
    n0 = n_elements(files)
    test = readfits(files[0],/silent)
    n1 = n_elements(test[*,0])
    n2 = n_elements(test[0,*])
    
    for j=0,n0-1 do begin
        data = readfits(files[j]+'.fits',/silent)
        mask = 'mask'+dat[j]
        ptow = 'ptow'+dat[j]
        readcol,silent=1,dotR[j],format='i,x,x,x,x,x',fiber
        bfi = where(fiber eq -1)
        readcol,silent=1,mask,format='i,x',y
        y = y-1
        for k=0,n_elements(bfi)-1 do $
          data[*,y[bfi[k]]-step:y[bfi[k]]+step,j] = -666
        outname = files[j]+'TC.fits'
        dataout = realign(datain,ptow,mask,wi,wdisp,step)
        writefits,outname,dataout
    endfor

endif    

END
