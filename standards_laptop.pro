; This code is used to extract lick standards from a collapsed data
; frame. This code works by entering a list of the following format:
; collapsed data frame
; ptow.dat file
; mask file
; the VIRUS-P_IR fits file for the month.
; the weight file
; HD number
; stellar type

; EX:
; jm3914_pefsmc.fits ptow_n1s.dat mask_may10_n1.dat fit_psfs_n1.fits jm3914_pefw.fits HD95272 K1_III

; The code requires a particular directory sub-structure. It will look
; for a sub-folder called "stars" where all the considerable output is
; dumped. It will also look for the flux_factor correction text in a
; specific directory.

; COMPILE strnsignif.pro
; .r strnsignif
; .r standards_laptop

PRO standards_laptop, list
;*********************************************************************
plotting = 'on'
threshold = 2000.0 ;The level to consider a fiber to contain a star
wi = 3550.0 
wdisp = 1.125
factor = 3.0 ;This is added to super-sample the data, in the same fashion as pipe2_v2
aperture = 5.0
step = (aperture - 1.0) / 2.0 
;slow = 'on' ;if set to 'on' a 2 sec delay occurs between plots
slow = 'off'
;*********************************************************************

readcol,list,f='a,a,a,a,a,a,a',files,ptowfiles,maskfiles,IRfiles,wgtfiles,star,type
n0 = n_elements(files)
test = readfits(files[0],/silent)
n1 = n_elements(test[*,0]) ; length of spectra in wavelength direction
n2 = n_elements(test[0,*])
waveO = dblarr(n1*factor)

for j=0,n1*factor-1 do waveO[j] = wi + ((wdisp/factor)*j)
i5000 = where(waveO ge 4950.0 and waveO lt 5050.0)

fmt0 = '(a10,2x,a8,2x,a6,2x,a9,2x,a10)'
ttl1 = ['   spectra','  km/sec','  FWHM','      S/N','wavelength']
fmt1 = '(f10.4,2x,f8.4,2x,f6.4,2x,f9.4,3x,f9.4)'

readcol,'/home/danu/murphy/research/flux_calib/fluxfactor_VP2_NFR_FI1_VPH2.txt',$
  silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
fluxcor = fluxcor * 1e-17

;added on sept 8th to smooth things down and see what we see...
sfluxcor = smooth(fluxcor,100,/edge)

window,2,retain=2
for j=0,n0-1 do begin
    frame = readfits(files[j],h,/silent)
    IR = readfits(IRfiles[j],/silent)
    wgt = readfits(wgtfiles[j],/silent)

    helioC = sxpar(h,'VHELIO',count=cnt)
    if (cnt eq 0) then helioC = 0.0
    print, 'The heliocentric correction is '+strn(helioC)
    helioC = 1.0 - (helioC / 299792.458) ;the variable 'helioC' is now a correction factor

    airmass = sxpar(h,'AIRMASS',cnt)
    if (cnt eq 0) then airmass = 1.0
    print,'The airmass is '+strn(airmass)
    pairmass = strnsignif(airmass,3)

    wavezp = sxpar(h,'WAVEZP',cnt)
    if (cnt eq 0) then wavezp = 0.0
    print, 'The wavelength zeropoint is '+strn(wavezp)

    readcol,ptowfiles[j],f='d,d,d,d,d',c0,c1,c2,c3,c4
    c0 = c0 + wavezp

    readcol,maskfiles[j],f='i,x',FL
    FLi = FL - 1

    fname = strsplit(files[j],'_.',/extract)
    fname = fname[0]

    chunk = frame[500:n1-200,*]
    med = median(chunk,dim=1)
    xa = findgen(n_elements(med)) + 1
    plot,xa,med,title=files[j],ytitle='CCD Counts',xtitle='Row',$
      xrange=[min(xa),max(xa)],/xstyle
    i = where(med gt threshold,count)
    print,'Number of fibers with a star in frame '+files[j]+': '+strn(count)
;    if (slow eq 'on') then wait,2 else wait,0.2
    
    if (count eq 0) then goto, jump1
  
    stars = med[i]
    istd = i[bsort(med[i])]
    istd = reverse(istd)
    Fflux = med[istd]
    Nflux = Fflux / max(Fflux)

    cd, 'stars'
    set_plot,'ps'
    device,file=fname+'_flux.ps'
    plot,xa,med,title=files[j]+'  '+star[j]+'  '+type[j],/xstyle,$
      ytitle='CCD Counts',xtitle='Fiber Number',xrange=[min(xa),max(xa)]
    device,/close_file
    set_plot,'x'
    cd, '../'
    if (slow eq 'on') then wait,2

    for k=0,count-1 do begin ;a loop through each found star spectra
        onespec = frame[*,istd[k]]
        oneIR = IR[*,istd[k]]
        onewgt = wgt[*,FLi[istd[k]]-step:FLi[istd[k]]+step]
        onewgt = total(onewgt,2) ;the weight array is collapsed
        oneflux = Nflux[k]
        wave = dblarr(n1)
        
        for l=0,n1-1 do wave[l] = c0[istd[k]] + c1[istd[k]]*l + c2[istd[k]]*l*l + $
          c3[istd[k]]*l*l*l + c4[istd[k]]*l*l*l*l
        wave = wave / helioC

        waveIR = IR[*,n_elements(IR[0,*])-1]

        i666 = where(onespec le -666.0,countbad)
        if (countbad gt 0) then onewgt[i666] = !Values.F_NAN
        if (countbad gt 0) then onespec[i666] = !Values.F_NAN

        A0 = interpol(onewgt,wave,waveO)  ;interpolated weights
        A1 = interpol(onespec,wave,waveO) ;interpolated spectra
        A2 = interpol(oneIR,waveIR,waveO)  ;interpolated instrumental dispersion (FWHM)

        s2n = fltarr(n1*factor)
        for l=0,n1*factor-1 do s2n[l] = sqrt(A0[l] * A1[l] * A1[l]) ;a small cheat here, as this works from collapsed data.

        A3 = s2n
        A4 = (A2 / waveO * 299792.458) / (2 * sqrt(2*alog(2))) ;Instrumental sigma

        s2n_est = round(mean(s2n[i5000]))
        corr1 = interpol(sfluxcor,fluxwave,waveO)
        corr2 = interpol(transcor,fluxwave,waveO)
        A1 = A1 * corr1 * (10^(0.4 * corr2 *airmass))
;        A1 = A1 / median(A1)
        A1 = A1 / mean(A1[i5000]) ;the spectra is normalized to the wavelength = 5000 A value
        if (k eq 0) then scale = A1 else A5 = A1 / scale

        namefits = star[j]+'_'+type[j]+'_s2n'+strn(s2n_est)+'_f'+strn(istd[k]+1)+'_X'+pairmass+'_'+fname+'.fits'
        nameps1 = star[j]+'_'+type[j]+'_s2n'+strn(s2n_est)+'_f'+strn(istd[k]+1)+'_X'+pairmass+'_'+fname+'.ps'
        nameps2 = star[j]+'_'+type[j]+'_s2n'+strn(s2n_est)+'_f'+strn(istd[k]+1)+'_X'+pairmass+'_'+fname+'N.ps'
        nametxt = star[j]+'_'+type[j]+'_s2n'+strn(s2n_est)+'_f'+strn(istd[k]+1)+'_X'+pairmass+'_'+fname+'.txt'

        if (plotting eq 'on') then begin
            if (k eq 0 and j eq 0) then window,0,retain=2
            if (k gt 0) then begin
                wset,0
                plot,waveO,A5,yrange=[0,3.0],xrange=[min(waveO),max(waveO)-100],/xstyle,/ystyle
            endif
            wset,2
            plot,waveO,A1,yrange=[0,2.0],position=[0.1,0.35,0.97,0.93],$
              title=star[j]+'   Stellar Type: '+type[j]+'  Fiber #'+strn(istd[k]+1)+'  Frame: '+fname,$
              xrange=[min(waveO),max(waveO)-100],/xstyle,ytitle='CCD counts',/ystyle,$
              xtickname=['','','','']
            plot,waveO,A2,position=[0.1,0.1,0.97,0.30],/noerase,yminor=4,$
              xrange=[min(waveO),max(waveO)-100],xtitle='Wavelength (A)',$
              ytitle='I.R. (FWHM)',/xstyle,yrange=[4,6],/ystyle,yticks=3
            xyouts,0.2,0.85,strn(oneflux),charsize=2.0,charthick=2,/normal

            if (slow eq 'on') then wait,1
            cd, 'stars'
            set_plot,'ps'
            device,file=nameps1
            plot,waveO,A1,yrange=[0,2.0],position=[0.14,0.35,0.97,0.93],$
              title=star[j]+'   Stellar Type: '+type[j]+'  Fiber #'+strn(istd[k]+1)+'  Frame: '+fname,$
              xrange=[min(waveO),max(waveO)-100],/xstyle,xthick=2,ythick=2,$
              xtickname=['','','',''],thick=2,charthick=2,ytitle='Normalized counts',/ystyle
            plot,waveO,A2,position=[0.14,0.1,0.97,0.30],/noerase,$
              xrange=[min(waveO),max(waveO)-100],xtitle='Wavelength (A)',$
              ytitle='I.R. (FWHM)',/xstyle,xthick=2,ythick=2,thick=2,$
              charthick=2,yrange=[4,6],/ystyle,yminor=4,yticks=3
            xyouts,0.2,0.85,strn(oneflux),charsize=1.5,charthick=2,/normal
            device,/close_file
   
            if (k gt 0) then begin
                device,file=nameps2
                plot,waveO,A5,yrange=[0,2.0],position=[0.14,0.35,0.97,0.93],$
                  title=star[j]+'   Stellar Type: '+type[j]+'  Fiber #'+strn(istd[k]+1)+'  Frame: '+fname,$
                  xrange=[min(waveO),max(waveO)-100],/xstyle,xthick=2,ythick=2,$
                  xtickname=['','','',''],thick=2,charthick=2,ytitle='Divided counts',/ystyle
                plot,waveO,A1,position=[0.14,0.1,0.97,0.30],/noerase,$
                  xrange=[min(waveO),max(waveO)-100],xtitle='Wavelength (A)',$
                  ytitle='Normalized Counts',/xstyle,xthick=2,ythick=2,thick=2,$
              charthick=2,yrange=[0,2],/ystyle,yminor=4,yticks=3
                xyouts,0.2,0.85,strn(oneflux),charsize=1.5,charthick=2,/normal
                device,/close_file
            endif
            
         set_plot,'x'
            cd, '../'
        endif

        axis1 = uint(n1 * factor)
        axis2 = uint(5)
        sxaddpar,headout,'SIMPLE','T'
        sxaddpar,headout,'BITPIX',-64,' Number of bits per data pixel'                  
        sxaddpar,headout,'NAXIS',2,' number of data axes'
        sxaddpar,headout,'NAXIS1',axis1,' old size * factor'
        sxaddpar,headout,'NAXIS2',axis2,''
        sxaddpar,headout,'EXTEND','T',' FITS data may contain extensions'
        sxaddpar,headout,'WaveINT',wi,' Initial wavelength'
        sxaddpar,headout,'WDISP',wdisp/factor,' Wavelength dispersion term'
        sxaddpar,headout,'S2N',s2n_est,' Estimate of the signal-2-noise (~5000A)'
        sxaddpar,headout,'AIRMASS',pairmass
        if (k eq 0 and j eq 0) then sxaddpar,headout,'COMMENT','ROW1: Wavelength'
        if (k eq 0 and j eq 0) then sxaddpar,headout,'COMMENT','ROW2: Signal-2-Noise'
        if (k eq 0 and j eq 0) then sxaddpar,headout,'COMMENT','ROW3: Instrumental Sigma (km/sec)'
        if (k eq 0 and j eq 0) then sxaddpar,headout,'COMMENT','ROW4: Instrumental Resolution (FWHM in A)'
        if (k eq 0 and j eq 0) then sxaddpar,headout,'COMMENT','ROW5: Normalized flux'
        
        out = [[waveO],[A3],[A4],[A2],[A1]];wavelength, s2n, Isigma, Ires, normalized flux 
        cd, 'stars'
        writefits,namefits,out,headout
        out = [[A1],[A4],[A2],[A3],[waveO]]
;        free_lun,5
;        openw,5,nametxt
;        printf,5,ttl1,format=fmt0
;        for m=0,n1*factor-1 do printf,5,out[m,0],out[m,1],out[m,2],out[m,3],out[m,4],format=fmt1
;        free_lun,5
        cd, '../'
    endfor
    jump1:
endfor
pause
while (!d.window ne -1) do wdelete,!d.window

stop
END
