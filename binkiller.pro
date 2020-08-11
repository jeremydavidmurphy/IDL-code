PRO binkiller, list
; This code is used to just kill all the extra bins created by
; bingen2D. It reads in a list of the wc b?n???? values. It then
; generates a cshell script called jj1 that you run to remove the
; extras.

readcol,list,f='i,x,x,a',wc,bins

n0 = n_elements(bins)
free_lun,10
openw,10,'jj1'

for j=0,n0-1 do begin
    wc1 = wc[j]
    sbin = strsplit(bins[j],'012',/extract)
    sbin = sbin[0]
    print,sbin
    if wc1 eq 1 then printf,10,'rm ',bins[j]
    if wc1 eq 2 then begin
        if sbin eq 'bin' then printf,10,'rm ',bins[j]
        if sbin eq 'bnn' then printf,10,'rm ',bins[j]
    endif
endfor
free_lun,10

stop
END
