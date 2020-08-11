pro rebin, file, binsize, type=type

; This code rebins spectra. It does nothing more than taking the mean,
; median or total of the given region.

if (n_elements(type) eq 0) then type = 'mean'

data = readfits(file)
n0 = n_elements(data)
print,'The number of elements in the spectra is '+strn(n0)

binned = fltarr(floor(n0/binsize))

if (n0 lt 32700) then begin

    cntr = 0
    pindex = 0
    for j=0,n0-binsize-1,binsize do begin
        if (type eq 'mean') then binned[cntr] = mean(data[j:j+binsize])
        if (type eq 'median') then binned[cntr] = median(data[j:j+binsize])
        if (type eq 'sum') then binned[cntr] = total(data[j:j+binsize])
        cntr = cntr+1
        pindex = [pindex,j]
    endfor
    pindex = pindex[1:*]
    binned = binned[0:n_elements(binned)-2]

endif else begin

    div = ceil(n0/32700.)
    cntr = 0
    binout = 0
    cntr2 = 0
    repeat begin
        
        piece = data[cntr2:cntr2+32699]
jump1:
        n1 = n_elements(piece)
        binned = fltarr(floor(n1/binsize))
        cntr1 = 0
        pindex = 0
        for j=0,n1-binsize-1,binsize do begin
            if (type eq 'mean') then binned[cntr1] = mean(piece[j:j+binsize])
            if (type eq 'median') then binned[cntr1] = median(piece[j:j+binsize])
            if (type eq 'sum') then binned[cntr1] = total(piece[j:j+binsize])
            cntr1 = cntr1+1
            pindex = double([pindex,j])
        endfor
        binned = binned[0:n_elements(binned)-2]
        binout = [binout,binned]
        pindex = pindex[1:*]
        if (cntr eq 0) then pinout = double(pindex) 
        if (cntr gt 0) and (cntr lt div) then begin
            pindex = pindex + max(pinout)
            pinout = double([pinout,pindex])
        endif
        cntr = cntr + 1
        cntr2 = cntr2 + 32700
        if (cntr eq div-1) then begin
            left = n0-cntr2-1
            piece = data[cntr2:cntr2+left]
            goto,jump1
        endif
    endrep until (cntr eq div)
    binned = binout[1:*]
    pindex = pinout[1:*]
    binned = binned[0:n_elements(binned)-2]
endelse

if (type ne 'total') then begin
    window,0,retain=2
    device,decomposed=0
    loadct,0
    plot,data
    loadct,4
    oplot,pindex,binned,color=150
endif

name = ''
print,'Name the output file:'
read,name

writefits,name,binned
wdelete
stop
end
