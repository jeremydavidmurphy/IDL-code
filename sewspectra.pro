PRO sewspectra, list

dataORs2n = 'data'
; This routine stitches together the spectra from pipe2_v2 into a
; single spectra.

; LIST: A list with two columns. The first column must be the blue
; spectra.

readcol,list,f='a,a',l1,l2
n0 = n_elements(l1)

for j=0,n0-1 do begin
    f1 = readfits(l1[j],h,/silent)
    f2 = readfits(l2[j],/silent)
    if (dataORs2n eq 'data') then begin
        lout = strsplit(l1[j],'_',/extract)
        lout = lout[0]+'.fits'
    endif else begin
        lout = strsplit(l1[j],'_',/extract)
        lout = lout[0]+'_S2N.fits'
    endelse
    out = [f1,f2]
    writefits,lout,out,h
endfor

stop
END
