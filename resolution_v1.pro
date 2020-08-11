; This code is used to estimate the instrumental resolution going into
; a given bin. It does this by fitting a gaussian to the arc lines for
; a given night and a given fiber. These values are collected into an
; array and then the median of all the FWHM are taken. This
; information is then used to convolve with the template stars for the
; kinematics.

PRO resolution, binlist, datalist

; BINLIST: This is a list of all the bins for a given data set. Each
; bin is a list of fibers and their pointing (123_b, 137_B, etc.)

; DATA.LIST: This is ALMOST the master list that goes into PIPE2.PRO. It is
; of the form
; FILE     POINTING     ?.DAT     dotR
; with just the pointing and it's corresponding ?.dat file that's read
; in. The difference is that it's been trimmed. This avoids multiple
; countings when several of the same exposures for a given pointing
; are taken on the same night.

; This routine calls the function rescheckF.pro, which is a function
; that does the gaussian fitting. It returns an n,m,2 array of the
; form:
; n: wavelength (11 of them- see below for their values)
; m: fiber number, with fiber #1 being at the top
; 2: FWHM and dispersion

;*****************************************************************
ForD = 0 ;use to return FWHM (ang) values
;ForD = 1 ;use to return instrumental dispersion (km/sec) values
;*****************************************************************

hgcd1 = [4046.55,4077.83,4358.33,4678.15,4799.91,4916.07,$
         5085.82,5154.66,5460.74]

readcol,binlist,format='a',blist
n0 = n_elements(blist) ; the number of bins

readcol,datalist,format='x,a,a,x',pointing,dat

ia = where(pointing eq 'a')
ib = where(pointing eq 'b')
ic = where(pointing eq 'c')
id = where(pointing eq 'd')
ie = where(pointing eq 'e')

for j=0,n0-1 do begin ; a loop through the bins
    onelist = blist[j] ; the name of the individual bin is selected
    print,'Working on bin '+onelist
    readcol,onelist,format='a',fiblist
    n1 = n_elements(fiblist) ;number of fibers in the bin list
    fibers = intarr(n1)
    letters = strarr(n1)
    for k=0,n1-1 do begin ;the bin list is parsed for pointing and fib # 
        temp = strsplit(fiblist[k],'_',/extract)
        fibers[k] = uint(temp[0])
        letters[k] = temp[1]
    endfor

    ifa = where(letters eq 'a')
    ifb = where(letters eq 'b')
    ifc = where(letters eq 'c')
    ifd = where(letters eq 'd')
    ife = where(letters eq 'e')

; here is a check to make sure things make sense and that you don't
; have pointings in your bin list and not in your data list...
    if (ia[0] eq -1 and ifa[0] ne -1) then stop
    if (ib[0] eq -1 and ifb[0] ne -1) then stop
    if (ic[0] eq -1 and ifc[0] ne -1) then stop
    if (id[0] eq -1 and ifd[0] ne -1) then stop
    if (ie[0] eq -1 and ife[0] ne -1) then stop

    if (ifa[0] ne -1) then begin
        nt = n_elements(ia)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                binaB = fltarr(3,n_elements(ifa))
                binaG = fltarr(3,n_elements(ifa))
                binaR = fltarr(3,n_elements(ifa))
            endif
            onearc = dat[ia[l]]
            out = rescheckF(onearc)
            binaB = [[binaB],[out[0:2,fibers[ifa]-1,ForD]]]
            binaG = [[binaG],[out[3:5,fibers[ifa]-1,ForD]]]
            binaR = [[binaR],[out[6:8,fibers[ifa]-1,ForD]]]
        endfor
        binaB = binaB[*,n_elements(ifa):*] ;the initial zeros are tossed
        binaG = binaG[*,n_elements(ifa):*]
        binaR = binaR[*,n_elements(ifa):*]
        na = n_elements(binaB[0,*])
    endif else na = 0
    
    if (ifb[0] ne -1) then begin
        nt = n_elements(ib)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'b' pointing involved
            if (l eq 0) then begin
                binbB = fltarr(3,n_elements(ifb))
                binbG = fltarr(3,n_elements(ifb))
                binbR = fltarr(3,n_elements(ifb))
            endif
            onearc = dat[ib[l]]
            out = rescheckF(onearc)
            binbB = [[binbB],[out[0:2,fibers[ifb]-1,ForD]]]
            binbG = [[binbG],[out[3:5,fibers[ifb]-1,ForD]]]
            binbR = [[binbR],[out[6:8,fibers[ifb]-1,ForD]]]
        endfor
        binbB = binbB[*,n_elements(ifb):*] ;the initial zeros are tossed
        binbG = binbG[*,n_elements(ifb):*]
        binbR = binbR[*,n_elements(ifb):*]
        nb = n_elements(binbB[0,*])
    endif else nb = 0
    
    if (ifc[0] ne -1) then begin
        nt = n_elements(ic)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'c' pointing involved
            if (l eq 0) then begin
                bincB = fltarr(3,n_elements(ifc))
                bincG = fltarr(3,n_elements(ifc))
                bincR = fltarr(3,n_elements(ifc))
            endif
            onearc = dat[ic[l]]
            out = rescheckF(onearc)
            bincB = [[bincB],[out[0:2,fibers[ifc]-1,ForD]]]
            bincG = [[bincG],[out[3:5,fibers[ifc]-1,ForD]]]
            bincR = [[bincR],[out[6:8,fibers[ifc]-1,ForD]]]
        endfor
        bincB = bincB[*,n_elements(ifc):*] ;the initial zeros are tossed
        bincG = bincG[*,n_elements(ifc):*]
        bincR = bincR[*,n_elements(ifc):*]
        nc = n_elements(bincB[0,*])
    endif else nc = 0
    
    if (ifd[0] ne -1) then begin
        nt = n_elements(id)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'd' pointing involved
            if (l eq 0) then begin
                bindB = fltarr(3,n_elements(ifd))
                bindG = fltarr(3,n_elements(ifd))
                bindR = fltarr(3,n_elements(ifd))
            endif
            onearc = dat[id[l]]
            out = rescheckF(onearc)
            bindB = [[bindB],[out[0:2,fibers[ifd]-1,ForD]]]
            bindG = [[bindG],[out[3:5,fibers[ifd]-1,ForD]]]
            bindR = [[bindR],[out[6:8,fibers[ifd]-1,ForD]]]
        endfor
        bindB = bindB[*,n_elements(ifd):*] ;the initial zeros are tossed
        bindG = bindG[*,n_elements(ifd):*]
        bindR = bindR[*,n_elements(ifd):*]
        nd = n_elements(bindB[0,*])
    endif else nd = 0
    
    if (ife[0] ne -1) then begin
        nt = n_elements(ie)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'e' pointing involved
            if (l eq 0) then begin
                bineB = fltarr(3,n_elements(ife))
                bineG = fltarr(3,n_elements(ife))
                bineR = fltarr(3,n_elements(ife))
            endif
            onearc = dat[ie[l]]
            out = rescheckF(onearc)
            bineB = [[bineB],[out[0:2,fibers[ife]-1,ForD]]]
            bineG = [[bineG],[out[3:5,fibers[ife]-1,ForD]]]
            bineR = [[bineR],[out[6:8,fibers[ife]-1,ForD]]]
        endfor
        bineB = bineB[*,n_elements(ife):*] ;the initial zeros are tossed
        bineG = bineG[*,n_elements(ife):*]
        bineR = bineR[*,n_elements(ife):*]
        nee = n_elements(bineB[0,*])
    endif else nee = 0
    nbin = na + nb + nc + nd + nee
    form1 = '(f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x)'
    formw = '(f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x)'
    free_lun,5
    openw,5,onelist+'.res'
    printf,5,hgcd1,format=formw
    if (na[0] ne 0) then for k=0,na-1 do $
      printf,5,[binaB[*,k],binaG[*,k],binaR[*,k]],format=form1
    if (nb[0] ne 0) then for k=0,nb-1 do $
      printf,5,[binbB[*,k],binbG[*,k],binbR[*,k]],format=form1
    if (nc[0] ne 0) then for k=0,nc-1 do $
      printf,5,[bincB[*,k],bincG[*,k],bincR[*,k]],format=form1
    if (nd[0] ne 0) then for k=0,nd-1 do $
      printf,5,[bindB[*,k],bindG[*,k],bindR[*,k]],format=form1
    if (nee[0] ne 0) then for k=0,nee-1 do $
      printf,5,[bineB[*,k],bineG[*,k],bineR[*,k]],format=form1
    free_lun,5
endfor

; The values are now written out into files named binname.res. These
; are all the FWHM (or I.RES) values for 9 wavelengths. The function
; resmedF.pro is called to trim these lists 

stop
end
