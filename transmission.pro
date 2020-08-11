PRO transmission, fiber

; This routine works off a set of two lists, each of which is a list
; of lists. The FIBER name given above gets superscripted in the
; following manner.
; If FIBER = 'a' then there must exist a list named a_all.list in the
; calling directory. This list is composed of a set of lists of the form:
;    a19_350.list
;    a19_365.list, etc.

;*****************************
; THIS IS CRITICAL (THE STRUCTURE OF THE NAME OF THE LIST) AS THIS IS
; WHERE THE WAVELENGTHS ARE PICKED OUT. IT CAN BE ANYTHING_WAVE.ANYTHING.
;
; Each of these lists is then of the form:
; Fib-a19-001-350.fit ; Fib-a19-001-350bg.fit
; Fib-a19-002-350.fit ; Fib-a19-002-350bg.fit
; Fib-a19-003-350.fit ; Fib-a19-003-350bg.fit
; OR
; Fib-a19-001-560.fit ; Fib-a19-001-560bg.fit
;                     ; Fib-a19-002-560bg.fit
;                     ; Fib-a19-003-560bg.fit

; and can be of any length. This allows for bad frames to be tossed,
; while the remainer are median combined.

; A MAXIMUM OF 8 FIBERS CAN BE REDUCED AT A GIVEN TIME

; METHOD: This can be set to 'mean', 'sum', or 'median'. It sets how
; the files will be combined before the calculations are made. Median
; should be used only when there are a minimum of three files FOR ALL
; DATA AND BACKGROUNDS. This is necessary as if the cosmics are median
; combined away in the data/background and not in the other, the
; results will be biased.
; NOTE: Using MEAN is very slow. ~30x longer than either sum or median

; All the baseline and data files need to exist in the same
; directory. 

;****************************************************************
; CHANGE THIS IF YOU'RE INTERESTED IN COMBINING THE DATA IN A
; DIFFERENT WAY...
METHOD = 'median'
;****************************************************************

readcol,fiber+'_all.list',format='a',fiblist
n0 = n_elements(fiblist)
wave = intarr(n0)

readcol,fiber+'_bg.list',format='a',baselist
n1 = n_elements(baselist)

if (n0 ne n1) then begin
    print,'You have different numbers of wavelengths in your BG and FIBER files!'
    stop
endif

for j=0,n0-1 do begin
    t = strsplit(fiblist[j],'_.',/extract)
    t = t[1]
    wave[j] = uint(t)
endfor

tbase = dblarr(n0,20,2)
tfib = dblarr(n0,20,2)
final = dblarr(n0,7)

values = dblarr(n_elements(wave),8)
out = dblarr(n_elements(wave),3,2)

window,1,retain=2,xsize=512,ysize=512,ypos=50
window,3,retain=2,xsize=512,ysize=512,ypos=50
loadct,37

for j=0,n0-1 do begin ; a loop through each wavelength
    flist = fiblist[j]
    blist = baselist[j]
    readcol,flist,format='a,a',delimiter=';',fibL,fibbgL
    readcol,blist,format='a,a',delimiter=';',baseL,basebgL

    nbi = where(baseL ne '')
    nb = n_elements(nbi)
    nbbi = where(basebgL ne '')
    nbb = n_elements(nbbi)
    nfi = where(fibL ne '')
    nf = n_elements(nfi)
    nfbi = where(fibbgL ne '')
    nfb = n_elements(nfbi)
    
    for k=0,n_elements(baseL)-1 do baseL[k]= strsplit(baseL[k],' ',/extract)
    for k=0,n_elements(basebgL)-1 do basebgL[k] = strsplit(basebgL[k],' ',/extract)
    for k=0,n_elements(fibL)-1 do fibL[k] = strsplit(fibL[k],' ',/extract)
    for k=0,n_elements(fibbgL)-1 do fibbgL[k] = strsplit(fibbgL[k],' ',/extract)

    base = dblarr(1024,1024,nb)
    basebg = dblarr(1024,1024,nbb)
    fib = dblarr(1024,1024,nf)
    fibbg = dblarr(1024,1024,nfb)

    for k=0,nb-1 do begin
        t = readfits(baseL[nbi[k]],silent=1)
        if (t[0] ne -1) then begin 
            base[*,*,k] = t
            tbase[j,k,0] = total(base[*,*,k])
            print,'The total for '+baseL[nbi[k]]+' is: '+strn(tbase[j,k,0])
        endif else begin
            print,'This file is in the list '+blist
            stop
        endelse
    endfor
    if (nb eq 1) then begin
        Mbase = base
        goto, jump1
    endif
    if (method eq 'median') then Mbase = median(base,dimension=3,/even)
    if (method eq 'mean') then begin
        Mbase = dblarr(1024,1024)
        for k=0,1023 do begin
            for l=0,1023 do begin
                Mbase[k,l] = mean(base[k,l,*])
            endfor
        endfor
    endif
    if (method eq 'sum') then Mbase = total(base,3)
    jump1:
    wset,3
    TVImage, BytScl(Mbase,Top=!D.Table_Size-3)
    
    for k=0,nbb-1 do begin
        t = readfits(basebgL[nbbi[k]],silent=1)
        if (t[0] ne -1) then begin
            basebg[*,*,k] = t
            tbase[j,k,1] = total(basebg[*,*,k])
            print,'The total for '+basebgL[nbbi[k]]+' is '+strn(tbase[j,k,1])
        endif else begin
            print,'This file is in the list '+blist
            stop
        endelse
    endfor
    if (nbb eq 1) then begin
        Mbasebg = basebg
        goto, jump2
    endif
    if (method eq 'median') then Mbasebg = median(basebg,dimension=3)
    if (method eq 'mean') then begin
        Mbasebg = dblarr(1024,1024)
        for k=0,1023 do begin
            for l=0,1023 do begin
                Mbasebg[k,l] = mean(basebg[k,l,*])
            endfor
        endfor
    endif
    if (method eq 'sum') then Mbasebg = total(basebg,3)
    jump2:
    
    for k=0,nf-1 do begin
        t = readfits(fibL[nfi[k]],silent=1)
        if (t[0] ne -1) then begin
            fib[*,*,k] = t
            tfib[j,k,0] = total(fib[*,*,k])
            print,'The total for '+fibL[[nfi[k]]]+' is: '+strn(tfib[j,k,0])
        endif else begin
            print,'This file is in the list '+flist
            stop
        endelse
    endfor
    if (nf eq 1) then begin
        Mfib = fib
        goto, jump3
    endif
    if (method eq 'median') then Mfib = median(fib,dimension=3)
    if (method eq 'mean') then begin
        Mfib = dblarr(1024,1024)
        for k=0,1023 do begin
            for l=0,1023 do begin
                Mfib[k,l] = mean(fib[k,l,*])
            endfor
        endfor
    endif
    if (method eq 'sum') then Mfib = total(fib,3)
    jump3:
    wset,1
    TVImage, BytScl(Mfib,Top=!D.Table_Size-3)    
    
    for k=0,nfb-1 do begin
        t = readfits(fibbgL[nfbi[k]],silent=1)
        if (t[0] ne -1) then begin
            fibbg[*,*,k] = t
            tfib[j,k,1] = total(fibbg[*,*,k])
            print,'The total for '+fibbgL[[nfbi[k]]]+' is '+strn(tfib[j,k,1])
        endif else begin
            print,'This file is in the list '+flist
            stop
        endelse
    endfor
    if (nfb eq 1) then begin
        Mfibbg = fibbg
        goto, jump4
    endif
    if (method eq 'median') then Mfibbg = median(fibbg,dimension=3)
    if (method eq 'mean') then begin
        Mfibbg = dblarr(1024,1024)
        for k=0,1023 do begin
            for l=0,1023 do begin
                Mfibbg[k,l] = mean(fibbg[k,l,*])
            endfor
        endfor
    endif
    if (method eq 'sum') then Mfibbg = total(fibbg,3)
    jump4:

    values[j,0] = total(Mbase)
    values[j,1] = total(Mbasebg)
    values[j,2] = values[j,0] - values[j,1] ;the subtracted baselines
    values[j,3] = total(Mfib)
    values[j,4] = total(Mfibbg)
    values[j,5] = values[j,3] - values[j,4]
    values[j,6] = values[j,3] / values[j,0]
    values[j,7] = values[j,5] / values[j,2] ;the subtracted, normalized values

    i = where(tbase[j,*,0] ne 0)
    tbase[j,19,0] = median(tbase[j,i,0],/even)
    i = where(tbase[j,*,1] ne 0)
    tbase[j,19,1] = median(tbase[j,i,1],/even)
    i = where(tfib[j,*,0] ne 0)
    tfib[j,19,0] = median(tfib[j,i,0],/even)
    i = where(tfib[j,*,1] ne 0)
    tfib[j,19,1] = median(tfib[j,i,1],/even)
    
    final[j,0] = tbase[j,19,0] ; the total median fiber values
    final[j,1] = tbase[j,19,1] ; the total median fiber background values
    final[j,2] = tfib[j,19,0]  ; the total median fiber values
    final[j,3] = tfib[j,19,1]  ; the total median fiber background values
    final[j,4] = (tfib[j,19,0]-tfib[j,19,1]) / (tbase[j,19,0]-tbase[j,19,1])
    final[j,5] = values[j,7] ; the values taken from images
    final[j,6] = values[j,6] ; the unsubtracted values taken from the images
endfor

wdelete,1,3
colors = [60,110,180,150,220,60,110,180,150,220,60,110,180,150,220,60,110,180,150,220,60,110,180,150,220]
symbols = [2,4,5,6,7,2,4,5,6,7,2,4,5,6,7,2,4,5,6,7,2,4,5,6,7,2,4,5,6,7,2,4,5,6,7]

;******************************************************************************************************
; PLOTS OF THE FIBER VALUES
;******************************************************************************************************

window,0,retain=2,xsize=420,ysize=300,xpos=425,ypos=375
device,decomposed=0
loadct,0
plot,wave,tfib[*,0,0],psym=1,title='Pre-subtracted Fiber Values',xrange=[min(wave)-20,max(wave)+20],xstyle=1,$
  /ynozero,symsize=2,xtitle='Wavelength (nm)',ytitle='Total Counts',charsize=1.2
loadct,4
cntr = 1
for j=0,8 do begin
    if (total(tfib[*,cntr,0]) ne 0) then begin
        oplot,wave,tfib[*,cntr,0],psym=symbols[cntr-1],color=colors[cntr-1],symsize=2
        cntr = cntr + 1
    endif
endfor
loadct,0
oplot,wave,tfib[*,19,0]

window,1,retain=2,xsize=420,ysize=300,xpos=425,ypos=55
device,decomposed=0
plot,wave,tfib[*,0,1],psym=1,title=' Fiber Background Values',xrange=[min(wave)-20,max(wave)+20],xstyle=1,$
  /ynozero,symsize=2,xtitle='Wavelength (nm)',ytitle='Total Counts',charsize=1.2
loadct,4
cntr = 1
for j=0,8 do begin
    if (total(tfib[*,cntr,1]) ne 0) then begin
        oplot,wave,tfib[*,cntr,1],psym=symbols[cntr-1],color=colors[cntr-1],symsize=2
        cntr = cntr + 1
    endif
endfor
loadct,0
oplot,wave,tfib[*,19,1]

;******************************************************************************************************
; PLOTS OF THE BASELINE VALUES
;******************************************************************************************************

window,3,retain=2,xsize=420,ysize=300,xpos=850,ypos=375
device,decomposed=0
plot,wave,tbase[*,0,0],psym=1,title='Pre-Subtracted Baseline Values',xrange=[min(wave)-20,max(wave)+20],xstyle=1,$
  /ynozero,symsize=2,xtitle='Wavelength (nm)',ytitle='Total Counts',charsize=1.2
loadct,4
cntr = 1
for j=0,8 do begin
    if (total(tbase[*,cntr,0]) ne 0) then begin
        oplot,wave,tbase[*,cntr,0],psym=symbols[cntr-1],color=colors[cntr-1],symsize=2
        cntr = cntr + 1
    endif
endfor
loadct,0
oplot,wave,tbase[*,19,0]

window,4,retain=2,xsize=420,ysize=300,xpos=850,ypos=55
device,decomposed=0
plot,wave,tbase[*,0,1],psym=1,title=' Baseline Background Values',xrange=[min(wave)-20,max(wave)+20],xstyle=1,$
  /ynozero,symsize=2,xtitle='Wavelength (nm)',ytitle='Total Counts',charsize=1.2
loadct,4
cntr = 1
for j=0,8 do begin
    if (total(tbase[*,cntr,1]) ne 0) then begin
        oplot,wave,tbase[*,cntr,1],psym=symbols[cntr-1],color=colors[cntr-1],symsize=2
        cntr = cntr +1
    endif
endfor

loadct,0
oplot,wave,tbase[*,19,1]
pause

window,0,retain=2,xsize=700,ysize=550,xpos=300,ypos=200
plot,wave,values[*,7],psym=2,xrange=[min(wave)-20,max(wave)+20],xstyle=1,$
  yrange=[0.0,1.0],ystyle=1,xtitle='Wavelength (nm)',ytitle='Normalized Transmission',$
  symsize=2,thick=2
oplot,wave,values[*,6],psym=4,symsize=1.0,thick=1.5 ;the unsubtracted values
loadct,4
oplot,wave,final[*,4],psym=1,color=150,symsize=2.0,thick=2.0

pause
wdelete,0,1,3,4

form1 = '(i13,1x,i13,1x,i13,1x,i13,1x,i13,1x,i13,1x,i13,1x,i13,1x,i13)'
form2 = '(f13.1,1x,f13.1,1x,f13.1,1x,f13.1,1x,f13.1,1x,f13.1,1x,f13.1,1x,f13.1,1x,f13.1)'
form3 = '(f13.10,1x,f13.10,1x,f13.10,1x,f13.10,1x,f13.10,1x,f13.10,1x,f13.10,1x,f13.10,1x,f13.10)'
openw,5,fiber+'_trans.txt'
printf,5,wave,format=form1
for j=0,3 do printf,5,final[*,j],format=form2
printf,5,final[*,6],format=form3
for j=4,5 do printf,5,final[*,j],format=form3
free_lun,5

set_plot,'x'
loadct,0
stop
end
