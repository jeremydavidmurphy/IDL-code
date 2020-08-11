; COMPILE: rescheckF
; COMPILE: biweightSF

; This code outputs the bin####.res files. They have the output in the
; form:
; wavelength   FWHM    weight
; The weight is the the sum of the number of fibers * number of
; exposures going into a given bin. For example. If there's 5
; exposures going into a frame and 7 fibers going into a bin then the
; weight is 5*7=35.

PRO resolution, binlist, datalist

; THIS CODE REQUIRES THAT THE DATA LIST BE FOR ONE MONTH/NIGHT
; ONLY. CURRENTLY THIS IS SET UP TO RETURN THE RESOLUTION FOR AN
; ENTIRE MONTH. IF YOU WANT TO GO NIGHT BY NIGHT THEN THE ONLY CHANGE
; TO MAKE SHOULD BE THE NAMING CONVENTION OF THE NAME OF THE DATALIST.

; The name of the datalist is important as it's used to name the
; output files. It needs to be of the form galaxy_MONTH/NIGHT.list.
; The code parses out the MONTH/NIGHT bit and writes it into the
; output file.

; BINLIST: This is a list of all the bins for a given data set. Each
; bin is a list of fibers and their pointing (123_b, 137_B, etc.)

; DATA.LIST: This is the master list that goes into PIPE1.PRO and
; PIPE2.PRO. It is of the form
; FILE     POINTING     ?.DAT     dotR
; with just the pointing and it's corresponding ?.dat file that's read
; in. This list DOES NOT GET TRIMMED. While multiple pointings get
; counted as such, this multiple counting acts as an impromtu weight
; in the calculation of the biweight, which sets the resolution
; value. The number of fibers and pointings then sets the weight that
; this piece of the instrumental resoultion will have. As different
; nights and months have different values, different resolution files
; are given weights which are taken into account when running rfitlov
; and rfitlovc.

; A switch below (ForD) switches from returning the resolution in
; angstroms or km/sec.

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

; Modified on Aug 8, 2009 to include the resmap routine.

pplot = 'yes' ;set to "yes" to plot the results for each bin and each night

hgcd1 = [4046.5539, 4077.8298, 4358.3253, 4678.1474,$
         4799.9080, 5085.8173, 5460.7366, 5769.6000, 5790.6239]

noffib = 246 ;change this if the number of fibers changes (i.e. with VP1)

; The following is a part of the ranges being selected for the final
; output files.
;**********************************************************************
; This section defines which arc lines are used in each wavelength
; region to estimate the resolution.
; modify this range and you need to change the other starred region!
wv1 = mean(hgcd1[0:1])
wv2 = mean(hgcd1[2:4])
wv3 = mean(hgcd1[5:6])
wv4 = mean(hgcd1[7:8])

;**********************************************************************

readcol,binlist,format='a',blist,silent=1
n0 = n_elements(blist) ; the number of bins

readcol,datalist,format='x,a,a,x',pointing,dat,silent=1

allarcs = fltarr(n_elements(hgcd1),noffib,n_elements(dat))

; the resolution is calculated for each arc in the data.list
for j=0,n_elements(dat)-1 do begin
    temp = rescheckF(dat[j])
    if (ForD eq 0) then temp = temp[*,*,0]
    if (ForD eq 1) then temp = temp[*,*,1]
    allarcs[*,*,j] = temp
endfor

ia1 = where(pointing eq 'a1')
ia2 = where(pointing eq 'a2')
ia3 = where(pointing eq 'a3')
ia  = where(pointing eq 'a')
ib  = where(pointing eq 'b')
ic  = where(pointing eq 'c')
id  = where(pointing eq 'd')
ie  = where(pointing eq 'e')
iif = where(pointing eq 'f')
ig  = where(pointing eq 'g')
ih  = where(pointing eq 'h')

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

    ifa1 = where(letters eq 'a1')
    ifa2 = where(letters eq 'a2')
    ifa3 = where(letters eq 'a3')
    ifa  = where(letters eq 'a')
    ifb  = where(letters eq 'b')
    ifc  = where(letters eq 'c')
    ifd  = where(letters eq 'd')
    ife  = where(letters eq 'e')
    iff  = where(letters eq 'f')
    ifg  = where(letters eq 'g')
    ifh  = where(letters eq 'h')

;POINTING A1 
    if (ia1[0] ne -1 and ifa1[0] ne -1) then begin
        nt = n_elements(ia1)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(ifa1))
            endif
            out = allarcs[*,fibers[ifa1]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(ifa1):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(ia1)*n_elements(ifa1)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_a1_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='A pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

;POINTING A2 
    if (ia2[0] ne -1 and ifa2[0] ne -1) then begin
        nt = n_elements(ia2)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(ifa2))
            endif
            out = allarcs[*,fibers[ifa2]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(ifa2):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(ia2)*n_elements(ifa2)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_a2_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='A2 pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

;POINTING A3 
    if (ia3[0] ne -1 and ifa3[0] ne -1) then begin
        nt = n_elements(ia3)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(ifa3))
            endif
            out = allarcs[*,fibers[ifa3]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(ifa3):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(ia3)*n_elements(ifa3)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_a3_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='A3 pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

;POINTING A 
    if (ia[0] ne -1 and ifa[0] ne -1) then begin
        nt = n_elements(ia)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(ifa))
            endif
            out = allarcs[*,fibers[ifa]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(ifa):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(ia)*n_elements(ifa)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_a_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='A pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

;POINTING B 
    if (ib[0] ne -1 and ifb[0] ne -1) then begin
        nt = n_elements(ib)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(ifb))
            endif
            out = allarcs[*,fibers[ifb]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(ifb):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(ib)*n_elements(ifb)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_b_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='B pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

;POINTING C
    if (ic[0] ne -1 and ifc[0] ne -1) then begin
        nt = n_elements(ic)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(ifc))
            endif
            out = allarcs[*,fibers[ifc]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(ifc):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(ic)*n_elements(ifc)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_c_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='C pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

;POINTING D 
    if (id[0] ne -1 and ifd[0] ne -1) then begin
        nt = n_elements(id)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(ifd))
            endif
            out = allarcs[*,fibers[ifd]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(ifd):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(id)*n_elements(ifd)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_d_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='D pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

;POINTING E 
    if (ie[0] ne -1 and ife[0] ne -1) then begin
        nt = n_elements(ie)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(ife))
            endif
            out = allarcs[*,fibers[ife]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(ife):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(ie)*n_elements(ife)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_e_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='E pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

;POINTING F 
    if (iff[0] ne -1 and iff[0] ne -1) then begin
        nt = n_elements(iff)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(iff))
            endif
            out = allarcs[*,fibers[iff]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(iff):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(iif)*n_elements(iff)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_f_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='F pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

;POINTING G 
    if (ig[0] ne -1 and ifg[0] ne -1) then begin
        nt = n_elements(ig)     ;the number of arcs going in
        for l=0,nt-1 do begin ;a loop through each 'a' pointing involved
            if (l eq 0) then begin
                bin = fltarr(n_elements(hgcd1),n_elements(ifg))
            endif
            out = allarcs[*,fibers[ifg]-1,ForD]
            bin = [[bin],[out]]
        endfor
        bin = bin[*,n_elements(ifg):*] ;the initial zeros are tossed
        na = n_elements(bin[0,*])
        bi1 = biweightSF(bin[0:1,*])
        bi2 = biweightSF(bin[2:4,*])
        bi3 = biweightSF(bin[5:6,*])
        bi4 = biweightSF(bin[7:8,*])
        wt = n_elements(ig)*n_elements(ifg)
        free_lun,5
        openw,5,onelist+'_'+month[1]+'_g_w'+strn(wt)+'.res'
        printf,5,wv1,bi1,wt,format = form2
        printf,5,wv2,bi2,wt,format = form2
        printf,5,wv3,bi3,wt,format = form2
        printf,5,wv4,bi4,wt,format = form2
        free_lun,5
        if (pplot eq 'yes') then begin
            if (na lt 95) then cstep = floor(195./na) else cstep = 1
            window,0,retain=2
            device,decomposed=0
            loadct,0
            plot,hgcd1,bin[*,0],yrange=[min(bin)-0.2,max(bin)+0.2],$
              ystyle=1,/nodata,title='G pointings for bin '+onelist,$
              xtitle='wavelength (A)',ytitle='FWHM'
            for k=0,na-1 do oplot,hgcd1,bin[*,k],psym=1
            loadct,4
            oplot,[wv1,wv2,wv3,wv4],[bi1,bi2,bi3,bi4],$
              psym=2,symsize=2,thick=2,color=150
            wait,1
        endif
    endif

endfor

stop
END
