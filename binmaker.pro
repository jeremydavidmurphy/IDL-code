;This routine is used to generate the final output spectra by summing
;over fibers within a given bin. ALL of the data for a given galaxy is
;read in and ALL the bins are generated at the same time. Both the
;data and the bins are read in as a list.

;*********************************************************************
pro binmaker, dlist, blist, SUFFIX=suf
;*********************************************************************

;The suffix is an addon to the output fits files. For example, if the
;name of the bin is:
;  bin1234
;Then the output will be:
;  bin1234.fits
;If the suffix keyword is used as:
;  SUFFIX='dark'
;Then the fits files will be output as:
;  bin1234_dark.fits

;The data must be in order (A,B,C,D, etc) in the list.
;EX: m87a.fits
;    m87b.fits
;    etc.

;The bin list is just that, a list of the bins.
;THIS CODE EXPECTS THAT EACH BIN LIST WILL CONTAIN A LIST OF THE
;FIBERS TO COMBINE IN THE FOLLOWING FORMAT
;EX:  124_a
;     137_a
;     123_b
;     etc. where the number is the fiber and the letter is the data
;     frame.

;THIS CODE CURRENTLY CAN HANDLE UP TO 6 POINTINGS FOR A GIVEN
;GALAXY.

readcol,dlist,format='A',files
readcol,blist,format='A',binlist

test = readfits(files[0])
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*])
n3 = n_elements(files)
n4 = n_elements(binlist)

data = dblarr(n1,n2,n3)
for l=0,n3-1 do data[*,*,l] = readfits(files[l])
;All the flags are set to zero so as not to confuse the sum.
data[where(data eq -666)] = 0.0

window,0,retain=2
loadct,0
for j=0,n4-1 do begin ;a loop thru all the bins
    readcol,binlist[j],format='A',bin
    n5 = n_elements(bin)
    sbin = strarr(2,n5)
;Next, the bin list is split to form a column of numbers and a column
;of letters.
    for k=0,n5-1 do sbin[*,k] = strsplit(bin[k],'_',/extract)
    ia = where(sbin[1,*] eq 'a') ;Index A
    ib = where(sbin[1,*] eq 'b')
    ic = where(sbin[1,*] eq 'c')
    id = where(sbin[1,*] eq 'd')
    ie = where(sbin[1,*] eq 'e')
    iff = where(sbin[1,*] eq 'f') ;IT WON'T LET ME CALL THIS "if"!
    nia = n_elements(ia) ;Number of Index for A, etc
    nib = n_elements(ib)
    nic = n_elements(ic)
    nid = n_elements(id)
    nie = n_elements(ie)
    nif = n_elements(iff)
    pa = dblarr(n1)
    pb = dblarr(n1)
    pc = dblarr(n1)
    pd = dblarr(n1)
    pe = dblarr(n1)
    pf = dblarr(n1)
    cntr = [0,0,0,0,0,0]
    
    if (nia ne 0) then begin
        if (nia gt 1) then begin
            fia = sbin[0,ia]-1 ;Fiber Index for the A's
            pa = data[*,fia,0]
            cntr[0] = 1
        endif else begin
            if (ia ne -1) then begin
                fia = sbin[0,ia]-1
                pa = data[*,fia,0]
                cntr[0] = 1
            endif
        endelse
    endif

    if (nib ne 0) then begin
        if (nib gt 1) then begin
            fib = sbin[0,ib]-1
            pb = data[*,fib,1]
            cntr[1] = 1
        endif else begin
            if (ib ne -1) then begin
                fib = sbin[0,ib]-1
                pb = data[*,fib,1]
                cntr[1] = 1
            endif
        endelse
    endif
    if (nic ne 0) then begin
        if (nic gt 1) then begin
            fic = sbin[0,ic]-1
            pc = data[*,fic,2]
            cntr[2] = 1
        endif else begin
            if (ic ne -1) then begin
                fic = sbin[0,ic]-1
                pc = data[*,fic,2]
                cntr[2] = 1
            endif
        endelse
    endif
    if (nid ne 0) then begin
        if (nid gt 1) then begin
            fid = sbin[0,id]-1
            pd = data[*,fid,3]
            cntr[3] = 1
        endif else begin
            if (id ne -1) then begin
                fid = sbin[0,id]-1
                pd = data[*,fid,3]
                cntr[3] = 1
            endif
        endelse
    endif
    if (nie ne 0) then begin
        if (nie gt 1) then begin
            fie = sbin[0,ie]-1
            pe = data[*,fie,3]
            cntr[4] = 1
        endif else begin
            if (ie ne -1) then begin
                fie = sbin[0,ie]-1
                pe = data[*,fie,3]
                cntr[4] = 1
            endif
        endelse
    endif
    if (nif ne 0) then begin
        if (nif gt 1) then begin
            fif = sbin[0,iff]-1
            pf = data[*,fif,3]
            cntr[5] = 1
        endif else begin
            if (iff ne -1) then begin
                fif = sbin[0,iff]-1
                pf = data[*,fif,3]
                cntr[5] = 1
            endif
        endelse
    endif
    
    out = [[pa],[pb],[pc],[pd],[pe],[pf]]
    
;--------------------------------------------------------------------
;The addition of a simple weighting based on the median of each fiber
;in the sum
;--------------------------------------------------------------------
    n6 = n_elements(out[0,*])
    wgt = dblarr(n6)
    for k=0,n6-1 do wgt[k] = median(out[*,k])
    wgt = wgt/max(wgt)
    
    for k=0,n6-1 do out[*,k] = out[*,k]*wgt[k]
;--------------------------------------------------------------------
    
    out = total(out,2)
    nsuf = n_elements(suf)
    if (nsuf eq 0) then begin
        name = binlist[j]+'.fits'
    endif else name = binlist[j]+'_'+suf+'.fits'
    
    writefits,name,out
    
    print,bin
    ymax = median(out[1500:1800,*])+100
    plot,out,title='Data for '+binlist[j],charsize=1.5,yrange=[0,ymax],$
      xrange=[0,1800],xstyle=1
pause
endfor

wdelete,0

STOP
end
