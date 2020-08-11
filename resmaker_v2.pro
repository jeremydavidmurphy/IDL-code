;obsolete- taken over by version 3
PRO resmaker_v2, datalist, binlist, fitsfile

; This routine uses the output of calcRES.pro to generate the
; bin####.res files needed to extract the kinematics via fitlov. This
; is an updated version of the original, which made some unpleasant
; simplifications. This is a modern version of resolution.pro

; The format needed for the rres routine is as follows. A list
; entitled res.list. Inside this is found the individual files that
; give the IR for a given bin. As a bin is composed of several fibers,
; each with their own IR values, resmaker_v2 outputs a ascii file of
; the following format:
; wavelength   FWHM    weight

; It is named bin####_f#_weight.res

; The weight is the the sum of the number of exposures going into a
; bin. For this reason resmaker_v2 reads in your datalist (the list
; needed by both pipe1 and pipe2). This allows the routine to know how
; many exposures go into a given fiber. It will also read in your
; bin.list.

; DATALIST: The name of the data list. i.e. the list you feed pipe1
; and pipe2.

; ex: IDL> resmaker_v2, 'M87f.list'

pslow = 'on' ;turn to 'off' to not slow the plotting routine down for visual inspection.
pntarr = ['a','b','c','d','e','f','g','a1','a2','a3'] ;the possible pointings. just add one to this string array, if needed
;**********************************************************************

readcol,binlist,f='a',blist,silent=1
n0 = n_elements(blist)

readcol,datalist,silent=1,f='x,a,a,x',pointings,dat
n1 = n_elements(dat)
np = n_elements(pntarr)

for j=0,n1-1 do begin
    t = strsplit(dat[j],'.',/extract)
    dat[j] = t[0]
endfor

for j=0,n0-1 do begin ;a loop through each bin.
    onebin = blist[j]
    print,'Working hard on bin '+onebin+'...'
    readcol,onebin,silent=1,f='a',temp,skip=1
    nfib = n_elements(temp)
    fibers = strarr(nfib) & ptng = fibers
    for k=0,nfib-1 do begin
        t = strsplit(temp[k],'_',/extract)
        fibers[k] = t[0]
        ptng[k] = t[1]
    endfor
    for k=0,np-1 do begin ;a loop over each element in the pntarr array (all the possible pointings)
        ipf = where(ptng eq pntarr[k],fibcnt) ;which fibers in the bin come from this pointing
        if (fibcnt gt 0) then begin ;the pointing exists in the bin
            ipp = where(pointings eq pntarr[k],ptcnt) ;which pointings in the datalist that get included
            ending = ''
            for l=0,ptcnt-1 do begin ;a loop through each of a specific pointings in the datalist
                starting = dat[l]
                if (starting eq ending) then begin
                    print,'Stepping over that redundancy!'
                    print,starting
                    goto,jumpover
                endif
                iin = where(dat eq starting,wgt) ;the index of the frames included
                IRfile = fitsfile
                oneframe = readfits(IRfile,header,/silent)
                n2 = n_elements(oneframe[0,*])
                n3 = n_elements(oneframe[*,0])
                if (n3 lt 20) then begin
                    wave = oneframe[*,n2-1]
                    oneframe = oneframe[*,1:*]
                endif else begin
                    wd = sxpar(header,'CDELT1')
                    wi = sxpar(header,'CRVAL1')
                    wave = findgen(n3)*wd + wi
                endelse
                ifiber = fibers[ipf]-1
                for fib=0,fibcnt-1 do begin ;a loop through each fiber
                    free_lun,5
                    openw,5,onebin+'_'+pointings[ipp[l]]+starting+'_f'+strn(fibers[fib])+'_w'+strn(wgt)+'.res'
                    slice = oneframe[*,ifiber[fib]]
                    n4 = n_elements(slice)
                    for p=0,n4-1 do printf,5,[strn(wave[p]),' ',strn(slice[p]),' ',strn(uint(wgt))]
                    free_lun,5
                endfor
                jumpover:
                ending = starting
            endfor
        endif ;if the pointing exists in the specific bin
    endfor ;a loop over each element in the pntarr array (i.e. the possible pointings)
endfor ;a loop over each bin

stop
END
