PRO resmaker_v3, datalist, binlist, fitsfile

; This routine uses the output of calcRES.pro to generate the
; bin####.res files needed to extract the kinematics via fitlov. This
; is an updated version of the original, which made some unpleasant
; simplifications. This is a modern version of resolution.pro

; The format needed for the rres routine is as follows. A list
; entitled res.list. Inside this is found the individual files that
; give the IR for a given bin. As a bin is composed of several fibers,
; each with their own IR values, resmaker_v3 outputs an ascii file of
; the following format:
; wavelength   FWHM    weight

; It is named bin####_f#_weight.res

; The weight is the the sum of the number of exposures going into a
; bin. For this reason resmaker_v3 reads in your datalist (the list
; needed by both pipe1 and pipe2). This allows the routine to know how
; many exposures go into a given fiber. It will also read in your
; bin.list.

; DATALIST: The name of the data list. i.e. the list you feed pipe1
; and pipe2.

; ex: IDL> resmaker_v2, 'M87f.list','bin.list','VIRUS-P_IR_FITALL_jun12.fits'

; MODIFIED ON FEB 27, 2012: This routine has been
; simplified. calcRES.pro now returns the IR for AN ENTIRE OBSERVING
; RUN. For this reason, this code no longer cares about pointings, but
; just fiber numbers. So, if a fiber # is used then it is assumed to
; have that same IR for all nights and all pointings it falls into. IF
; YOU HAVE DATA THAT YOU NEED TO COMBINE FROM DIFFERENT OBSERVING RUNS
; AND YOU FEEL THE IR HAS CHANGED, YOU NEED TO RUN THIS CODE ON
; INDIVIDUAL OBSERVING RUNS, THEN APPLY ALL THE RESULTING BIN#.RES
; FILES DURING THE KINEMATIC EXTRACTION.

;pntarr = ['a','b','c','d','e','f','g','a1','a2','a3'] ;the possible pointings. just add one to this string array, if needed
;**********************************************************************

readcol,binlist,f='a',blist,silent=1
n0 = n_elements(blist)

readcol,datalist,silent=1,f='x,a,x,x,x,x',pointings

month = strsplit(fitsfile,'_.',/extract)
month = month[3]

IRfile = fitsfile
IRframe = readfits(IRfile,header,/silent)
n2 = n_elements(IRframe[0,*])
n3 = n_elements(IRframe[*,0])

if (n3 lt 20) then begin
    wave = IRframe[*,n2-1]
    IRframe = IRframe[*,1:*]
endif else begin
    wd = sxpar(header,'CDELT1')
    wi = sxpar(header,'CRVAL1')
    wave = findgen(n3)*wd + wi
endelse

for j=0,n0-1 do begin ;a loop through each bin.
    onebin = blist[j]
    print,'Working hard on bin '+onebin+'...'
    readcol,onebin,silent=1,f='a',temp
    i = where(temp ne -1)
    temp = temp[i]
    nfib = n_elements(temp)
    fibers = strarr(nfib) & ptng = fibers
    donefibers = 0
    weight = 0
    for k=0,nfib-1 do begin
        t = strsplit(temp[k],'_',/extract)
        fibers[k] = t[0]
        ptng[k] = t[1]
    endfor
    for k=0,nfib-1 do begin ;a loop over each fiber in the bin list
        weight = 0
        onefib = fibers[k]
        ifib = onefib-1
        i = where(onefib eq donefibers,cnt)
        if (cnt gt 0) then goto,jumpfiber ;if a fiber has been taken care of
        ii = where(onefib eq fibers,nf) ;determines whether the same fiber is used in multiple pointings (i.e. for a dither set)
        for l=0,nf-1 do begin ;loop over all the times the fiber is used for a different pointing
            onept = ptng[ii[l]]
            i = where(onept eq pointings,npt) ;determines the number of times the pointing for a fiber exists in the data.list
            if (npt lt 1) then begin
                print,'You have a pointing in your bin that is not in your data list!'
                print,onebin
                stop
            endif
            weight = weight + npt
        endfor
        
        free_lun,5 ;the results are written out
        openw,5,onebin+'_'+month+'_f'+strn(onefib)+'_w'+strn(weight)+'.res'
        slice = IRframe[*,ifib]
        n4 = n_elements(slice)
        for p=0,n4-1 do printf,5,[strn(wave[p]),' ',strn(slice[p]),' ',strn(uint(weight))]
        free_lun,5

        jumpfiber:
        donefibers = [donefibers,onefib]
    endfor
endfor ;a loop over each bin

stop
END
