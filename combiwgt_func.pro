;The key piece of pipe2. This makes a biweight combination of all the
;frames for a given pointing. The dotR files are the output of pipe1
;and serve to mask the foreground stars, etc, that get into a
;frame. If more than 1/2 of the frames to be combined have a specific
;fiber that's been masked, then the assumption is made that there is
;an object and the entire fiber get's a -666 flag.

;This routine then calls the function biweight.pro, after the proper
;fibers have been masked.

function combiwgt_func, array, dotR, NAMES=names,WEIGHT=wgt

;If the weight array is NOT used then the weights are calculated by
;taking the median of each data frame and each fiber, then normalizing
;by the median value of all data files for a given fiber.

;The NAME keyword can be used to write the .fits files out AFTER they
;have been masked. These files amount to the actual files used in the
;combination. They are named "name_AM.fits", A for aligned and M for
;masked. If the realign keyword was set to FREE, then they have also
;been trimmed to be of the same length for all frames, unlike the
;name_A.fits files written previously in PIPE2.PRO which may vary,
;depending upon the wavelength solution found in the ptow.dat file.

;The dotR is A LIST OF THE NAME.R FILES, not the actual files.
;The dotR arrays are of the form:
;fiber #    radius     deltaRA    deltaDec    normalized flux
;Yet it's only whether the fiber exists in the first column or not
;that matters.

n1 = n_elements(array[*,0,0])
n2 = n_elements(array[0,*,0])
n3 = n_elements(array[0,0,*])

if (n_elements(wgt) eq 0) then begin
    swch = 'off'
    wgt = dblarr(n3)
    wgtarray = dblarr(n3,n2)
    print,'The weights will calculated FOR EACH FIBER INDIVIDUALLY!'
endif else swch = 'on'

biwtout = dblarr(n1,n2)
fibarr = intarr(n3,n2)
slice = dblarr(n1,n3)

;The fiber rejection (based on the dotR files) is begun.
for j=0,n3-1 do begin
    readcol,dotR[j],format='I,X,X,X,X',gf
    gfi = gf[bsort(gf)]-1
    af = indgen(n2)+1
    af[gfi]=-1 ;here the GOOD fibers get a -1 tag
    fibarr[j,*] = af
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;This loop solves(?) the issue of flagging a majority, but not all
;pointings for a given fiber. If the majority of the exposures for a
;fiber are good then it stands, if over 1/2 are rejected, then the
;entire fiber is rejected. This should handle the low-lying stars that
;get noticed in only some of the exposures.
for k=0,n2-1 do begin
    m = median(fibarr[*,k])
    if (m ne -1) then fibarr[*,k] = 0
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;The rejected fibers (as per the name.R file) are tossed from the data
for l=0,n3-1 do begin
    for k=0,n2-1 do begin
        if (fibarr[l,k] ne -1) then begin
            array[*,k,l] = -666
;            print,'Fiber '+strn(k+1)+' in frame '+strn(l+1)+' was rejected!'
        endif
    endfor
endfor

if (n_elements(names) ne -1) then begin
    for l=0,n3-1 do begin
        name = names[l]+'_AM.fits'
        writefits,name,array[*,*,l]
    endfor
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for k=0,n2-1 do begin           ;This loops on each fiber
    print,'Working on fiber '+strn(k+1)
;all the exposures for one fiber are defined
    for l=0,n3-1 do slice[*,l] = array[*,k,l]
    if (swch eq 'off') then begin
         for l=0,n3-1 do begin
            wgt[l] = median(slice[*,l])
            if (wgt[l] lt 0.0) then wgt[l] = 0.0
        endfor
        t = total(wgt)
        if (t ne 0) then begin
            wgt = wgt/max(wgt)
            wgtarray[*,k] = wgt
        endif else begin
            biwtout[*,k] = -666
            wgtarray[*,k] = intarr(n3)
            goto,jumpend
        endelse
    endif

    biwtout[*,k] = biweight(slice,NORM=wgt)

jumpend:
endfor

writefits,'weights.fits',wgtarray
return,biwtout
end

