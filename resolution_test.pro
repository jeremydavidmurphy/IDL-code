; The change to this version from resolution_v1 is that now the
; datalist is for a given night rather than all the data for all the
; nights. This allows the resolution for a given (night?/run?) be
; calculated separately. 

PRO resolution, binlist, datalist

; THIS CODE REQUIRES THAT THE DATA LIST BE FOR ONE MONTH/NIGHT
; ONLY. CURRENTLY THIS IS SET UP TO RETURN THE RESOLUTION FOR AN
; ENTIRE MONTH. IF YOU WANT TO GO NIGHT BY NIGHT THEN THE ONLY CHANGE
; TO MAKE SHOULD BE THE NAMING CONVENTION. 

; BINLIST: This is a list of all the bins for a given data set. Each
; bin is a list of fibers and their pointing (123_b, 137_B, etc.)

; DATA.LIST: This is ALMOST the master list that goes into PIPE2.PRO. It is
; of the form
; FILE     POINTING     ?.DAT     dotR
; with just the pointing and it's corresponding ?.dat file that's read
; in. The difference is that it's been trimmed. This avoids multiple
; countings when several of the same exposures for a given pointing
; are taken on the same night.

; THIS MUST BE OF THE FORM NAME_MONTH_WHATEVER...

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

pplot = 'yes' ;set to "yes" to plot the results as they come

hgcd1 = [4046.55, 4077.83, 4358.33, 4678.15, 4799.91, 4916.07,$
         5085.82, 5460.74, 5769.60, 5790.62]
hgcdf = [4046.5539, 4077.8298, 4358.3253, 4678.1474, 4799.9080, 4916.0661,$
         5085.8173, 5460.7366, 5769.6000, 5790.6239]
; The following is a part of the ranges being selected for the final
; output files.
;**********************************************************************
; modify this range and you need to change the other starred region!
wv1 = mean(hgcdf[0:1])
wv2 = mean(hgcdf[2:4])
wv3 = mean(hgcdf[5:6])
wv4 = mean(hgcdf[7:9])
;**********************************************************************

readcol,binlist,format='a',blist,silent=1
n0 = n_elements(blist) ; the number of bins

readcol,datalist,format='x,a,a,x',pointing,dat,silent=1

ia = where(pointing eq 'a')
ib = where(pointing eq 'b')
ic = where(pointing eq 'c')
id = where(pointing eq 'd')
ie = where(pointing eq 'e')

month = strsplit(datalist,'_.',/extract)
form2 = '(f8.2,1x,f8.4,1x,i4)'

for j=0,n0-1 do begin ; a loop through the bins
    onelist = blist[j] ; the name of the individual bin is selected
    print,''
    print,'Working on bin '+onelist
    readcol,onelist,format='a',fiblist,silent=1
    print,'The fibers are...'
    print,fiblist
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

    biA = 0.0 & biB = 0.0 & biC = 0.0 & biD = 0.0 & biE = 0.0

    if (ia[0] ne -1 and ifa[0] ne -1) then begin
        nt = n_elements(ia)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                binA = fltarr(n_elements(hgcd1),n_elements(ifa))
            endif
            onearc = dat[ia[l]]
            out = rescheckF(onearc)
            binA = [[binA],[out[*,fibers[ifa]-1,ForD]]]
        endfor
        binA = binA[*,n_elements(ifa):*] ;the initial zeros are tossed
        na = n_elements(binA[0,*])
        biA1 = biweightSF(binA[0:1,*])
        biA2 = biweightSF(binA[2:4,*])
        biA3 = biweightSF(binA[5:6,*])
        biA4 = biweightSF(binA[7:9,*])
        wtA = n_elements(ia)*n_elements(ifa)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_w'+strn(wtA)+'.res'
        printf,5,wv1,biA1,format = form2
        printf,5,wv2,biA2,format = form2
        printf,5,wv3,biA3,format = form2
        printf,5,wv4,biA4,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,binA[*,0],yrange=[min(binA)-0.2,max(binA)+0.2],$
              ystyle=1,/nodata,title='A pointings for bin '+onelist,$
              xtitle='FWHM',ytitle='Wavelength (A)'
            for k=0,na-1 do oplot,hgcd1,binA[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[biA1,biA2,biA3,biA4],psym=2,symsize=2,thick=2,color=150
            pause
        endif
    endif else na = 0
    
    if (ib[0] ne -1 and ifb[0] ne -1) then begin
        nt = n_elements(ib)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'b' pointing involved
            if (l eq 0) then begin
                binB = fltarr(n_elements(hgcd1),n_elements(ifb))
            endif
            onearc = dat[ib[l]]
            out = rescheckF(onearc)
            binB = [[binB],[out[*,fibers[ifb]-1,ForD]]]
        endfor
        binB = binB[*,n_elements(ifb):*] ;the initial zeros are tossed
        nb = n_elements(binB[0,*])
        biB1 = biweightSF(binB[0:1,*])
        biB2 = biweightSF(binB[2:4,*])
        biB3 = biweightSF(binB[5:6,*])
        biB4 = biweightSF(binB[7:9,*])
        wtB = n_elements(ib)*n_elements(ifb)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_w'+strn(wtB)+'.res'
        printf,5,wv1,biB1,format = form2
        printf,5,wv2,biB2,format = form2
        printf,5,wv3,biB3,format = form2
        printf,5,wv4,biB4,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            window,1,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,binB[*,0],yrange=[min(binB)-0.2,max(binB)+0.2],$
              ystyle=1,/nodata,title='B pointings for bin '+onelist,$
              xtitle='FWHM',ytitle='Wavelength (A)'
            for k=0,nb-1 do oplot,hgcd1,binB[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[biB1,biB2,biB3,biB4],psym=2,symsize=2,thick=2,color=150
        endif
    endif else nb = 0

    if (ic[0] ne -1 and ifc[0] ne -1) then begin
        nt = n_elements(ic)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'c' pointing involved
            if (l eq 0) then begin
                binC = fltarr(n_elements(hgcd1),n_elements(ifc))
            endif
            onearc = dat[ic[l]]
            out = rescheckF(onearc)
            binC = [[binC],[out[*,fibers[ifc]-1,ForD]]]
        endfor
        binC = binC[*,n_elements(ifc):*] ;the initial zeros are tossed
        nc = n_elements(binC[0,*])
        biC1 = biweightSF(binC[0:1,*])
        biC2 = biweightSF(binC[2:4,*])
        biC3 = biweightSF(binC[5:6,*])
        biC4 = biweightSF(binC[7:9,*])
        wtC = n_elements(ic)*n_elements(ifc)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_w'+strn(wtC)+'.res'
        printf,5,wv1,biC1,format = form2
        printf,5,wv2,biC2,format = form2
        printf,5,wv3,biC3,format = form2
        printf,5,wv4,biC4,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            window,2,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,binC[*,0],yrange=[min(binC)-0.2,max(binC)+0.2],$
              ystyle=1,/nodata,title='C pointings for bin '+onelist,$
              xtitle='FWHM',ytitle='Wavelength (A)'
            for k=0,nc-1 do oplot,hgcd1,binC[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[biC1,biC2,biC3,biC4],psym=2,symsize=2,thick=2,color=150
        endif
    endif else nc = 0
    
    if (id[0] ne -1 and ifd[0] ne -1) then begin
        nt = n_elements(id)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'd' pointing involved
            if (l eq 0) then begin
                binD = fltarr(n_elements(hgcd1),n_elements(ifd))
            endif
            onearc = dat[id[l]]
            out = rescheckF(onearc)
            binD = [[binD],[out[*,fibers[ifd]-1,ForD]]]
        endfor
        binD = binD[*,n_elements(ifd):*] ;the initial zeros are tossed
        nd = n_elements(binD[0,*])
        biD1 = biweightSF(binD[0:1,*])
        biD2 = biweightSF(binD[2:4,*])
        biD3 = biweightSF(binD[5:6,*])
        biD4 = biweightSF(binD[7:9,*])
        wtD = n_elements(id)*n_elements(ifd)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_w'+strn(wtD)+'.res'
        printf,5,wv1,biD1,format = form2
        printf,5,wv2,biD2,format = form2
        printf,5,wv3,biD3,format = form2
        printf,5,wv4,biD4,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            window,3,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,binD[*,0],yrange=[min(binD)-0.2,max(binD)+0.2],$
              ystyle=1,/nodata,title='D pointings for bin '+onelist,$
              xtitle='FWHM',ytitle='Wavelength (A)'
            for k=0,nd-1 do oplot,hgcd1,binD[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[biD1,biD2,biD3,biD4],psym=2,symsize=2,thick=2,color=150
        endif
    endif else nd = 0
    
    if (ie[0] ne -1 and ife[0] ne -1) then begin
        nt = n_elements(ie)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'e' pointing involved
            if (l eq 0) then begin
                binE = fltarr(n_elements(hgcd1),n_elements(ife))
            endif
            onearc = dat[ie[l]]
            out = rescheckF(onearc)
            binE = [[binE],[out[*,fibers[ife]-1,ForD]]]
        endfor
        binE = binE[*,n_elements(ife):*] ;the initial zeros are tossed
        nee = n_elements(binE[0,*])
        biE1 = biweightSF(binE[0:1,*])
        biE2 = biweightSF(binE[2:4,*])
        biE3 = biweightSF(binE[5:6,*])
        biE4 = biweightSF(binE[7:9,*])
        wtE = n_elements(ie)*n_elements(ife)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_w'+strn(wtE)+'.res'
        printf,5,wv1,biE1,format = form2
        printf,5,wv2,biE2,format = form2
        printf,5,wv3,biE3,format = form2
        printf,5,wv4,biE4,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            window,4,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,binE[*,0],yrange=[min(binE)-0.2,max(binE)+0.2],$
              ystyle=1,/nodata,title='E pointings for bin '+onelist,$
              xtitle='FWHM',ytitle='Wavelength (A)'
            for k=0,nee-1 do oplot,hgcd1,binE[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[biE1,biE2,biE3,biE4],psym=2,symsize=2,thick=2,color=150
        endif
    endif else nee = 0
endfor


stop
END
