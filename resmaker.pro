PRO resmaker, IRfits

; This routine is taking over the resolution.pro routine. The output
; of calcRES.pro is now a fits file of the IR for a month. If all the
; data for a galaxy is from a single month, then 
; you've got a single IR map. This is a simpler code to generate the
; bin.res files needed for the Monte Carlo simulations.

;**********************************************************************
MorMorA = 'all' ;set this to 'median' or 'max' depending on whether you want to return the median of the IR or and max for each given spatial bin.
; set to 'all' if you want to output the wavelength, min, median, max
pslow = 'on' ;turn to 'off' to not slow the plotting routine down for visual inspection.
;**********************************************************************

readcol,'bin.list',f='a',bins,silent=1
n0 = n_elements(bins)

IRarr = readfits(IRfits,/silent)
n1 = n_elements(IRarr[0,*])

wave = IRarr[*,n1-1]
IRarr = IRarr[*,0:n1-2]
n2 = n_elements(wave)

window,0,retain=2

for j=0,n0-1 do begin
    bin = bins[j]
    print,'Working on bin ' +bin
    readcol,bin,silent=1,f='a',fibers
    ikeep = where(fibers ne -1,count)
    if (count gt 0) then fibers = fibers[ikeep]
    n3 = n_elements(fibers)
    ifiber = intarr(n3)
    for k=0,n3-1 do begin
        a = strsplit(fibers[k],'_',/extract)
        ifiber[k] = uint(a[0]) - 1.0
    endfor
    IRone = fltarr(n2,n3)
    for k=0,n3-1 do IRone[*,k] = IRarr[*,ifiber[k]]
    if (MorMorA eq 'median') then IRbin = median(IRone,dim=2,/even)
    if (MorMorA eq 'max') then begin
        IRbin = fltarr(n2)
        for k=0,n2-1 do IRbin[k] = max(IRone[k,*])
    endif
    if (MorMorA eq 'all') then begin
        IRbin = fltarr(n2,3)
        for k=0,n2-1 do begin
            IRbin[k,0] = min(IRone[k,*])
            IRbin[k,1] = median(IRone[k,*],/even)
            IRbin[k,2] = max(IRone[k,*])
        endfor
    endif
    plot,wave,IRone[*,0],/nodata,yrange=[min(IRone),max(IRone)],$
      title=bin,xtitle='Wavelength (A)',ytitle='FWHM (A)'
    for k=0,n3-1 do oplot,wave,IRone[*,k],psym=1
    oplot,wave,IRbin,thick=2
    if (bin eq 'bin1505') then begin
        set_plot,'ps'
        device,file='bin1505_res.ps'
        plot,wave,IRone[*,0],/nodata,yrange=[min(IRone),max(IRone)],$
          title=bin,xtitle='Wavelength (A)',ytitle='FWHM (A)',$
          xthick=2,ythick=2,charthick=2,thick=2
        for k=0,n3-1 do oplot,wave,IRone[*,k],psym=1,thick=2
        device,/close_file
        set_plot,'x'
    endif
        
    if (pslow eq 'on') then wait,0.5
    
    free_lun,5
    openw,5,bin+'.res'
    if (MorMorA ne 'all') then begin
        for k=0,n2-1 do printf,5,[strn(wave[k]),'  ',strn(IRbin[k]),'  100']
        free_lun,5
    endif
    if (MorMorA eq 'all') then begin
        for k=0,n2-1 do printf,5,[strn(wave[k]),'  ',strn(IRbin[k,0]),'  ',$
                                  strn(IRbin[k,1]),'  ',strn(IRbin[k,2])]
        free_lun,5
    endif
endfor

stop
END
