PRO binresmaker, datalist, binlist, fitslist

; This code is used to determine the median and S.D. of the
; instrumental resolution (IR) going into a given spatial bin. The
; code uses the output of calcRES.pro (the FITALL fits files), the bin
; list (bin.list) and the data list (GALNAMEdata.list) to generate an
; output files file that has the wavelength, IR and IR uncertainty.

; DATALIST: Just your data list (e.g. M87data.list) from pipe2

; BINLIST: Just your bin.list

; FITSLIST: A list of two columns, the observing run (e.g. jan08,
; march12, etc) and the VIRUS-P IR fits file (e.g. VIRUS-P_IR_FITALL_jan08.fits).

ploton = 'on'
;**********************************************************************

readcol,binlist,f='a',blist,silent=1
n0 = n_elements(blist)

;readcol,datalist,silent=1,f='x,a,a,x,x,x,x,x',pointings,runs
readcol,datalist,silent=1,f='x,a,a,x,x,x',pointings,runs ;the the older lists, before the sky noise was part of the routine...

n1 = n_elements(pointings)

readcol,fitslist,silent=1,f='a,a',obsruns, irfiles
n2 = n_elements(irfiles)

temp = readfits(irfiles[0],/silent,h)
n3 = n_elements(temp[*,0])
wd = sxpar(h,'CDELT1')
wi = sxpar(h,'CRVAL1')
wave = findgen(n3)*wd + wi

sxaddpar,h,'COMMENT','Row 1: Mean Instrumental Resolution (Ang, FWHM)'
sxaddpar,h,'COMMENT','Row 2: SD Instrumental Resolution (Ang, FWHM)'
sxaddpar,h,'COMMENT','Row 3: Min Instrumental Resolution (Ang, FWHM)'
sxaddpar,h,'COMMENT','Row 4: Max Instrumental Resolution (Ang, FWHM)'
sxaddpar,h,'COMMENT','Row 5: Wavelength (A)'
countplot = 0
ans = ''

for j=0,n0-1 do begin ;a loop through each bin.
   irarray = fltarr(n3)
    onebin = blist[j]
    print,'Working hard on bin '+onebin+'...'
    readcol,onebin,silent=1,f='a',temp
    i = where(temp ne -1)
    temp = temp[i]
    nfib = n_elements(temp)
    fibers = strarr(nfib) & ptng = fibers
    donefibers = 0
    weight = 0
    for k=0,nfib-1 do begin
        t = strsplit(temp[k],'_',/extract)
        fibers[k] = t[0]
        ptng[k] = t[1]
     endfor
    ifib = fibers - 1
    for k=0,nfib-1 do begin ;a loop over each fiber in the bin list
       print,'Working on fiber '+strn(k+1)+'...'
       ipoint = where(pointings eq ptng[k],countpoint)
       if (countpoint eq 0) then begin
          print,''
          print,'That pointing is not in your data list!!!'
          goto,jumpfiber
       endif
       maskfiles = runs[ipoint]
       for l=0,countpoint-1 do begin ;a loop through the number of exposures
          temp = strsplit(runs[ipoint[l]],'_',/extract)
          irun = where(obsruns eq temp[0],countrun)
          if (countrun eq 0) then begin
             print,''
             print,'The IR map was not found for one of your pointings!!!'
             stop
          endif
          ir = readfits(irfiles[irun],/silent)
       irarray = [[irarray],[ir[*,ifib[k]]]]
       endfor ;end of the loop for each night
       jumpfiber: 
    endfor ;end of the fiber loop
    irarray = irarray[*,1:*]
    ibad = where(irarray le 0.0,countbad)
    if (countbad gt 0) then irarray[ibad] = !VALUES.F_NAN

    irout = fltarr(n3,5)
    irout[*,4] = wave
    for w=0,n_elements(irarray[*,0])-1 do begin
       irout[w,0] = mean(irarray[w,*],/nan)
       irout[w,1] = stddev(irarray[w,*],/nan)
       irout[w,2] = min(irarray[w,*])
       irout[w,3] = max(irarray[w,*])
    endfor

    writefits,onebin+'_IR.fits',irout,h

    if (ploton eq 'on') then begin
       if (countplot eq 0) then begin
          window,0,retain=2
          countplot = 1
       endif
       loadct,0,/silent
       plot,irout[*,4],irout[*,0]+irout[*,1],linestyle=3,ytitle='Instrumental Resolution',$
            xtitle='Wavelength (A)',title=onebin,/nodata,yrange=[4,7],ys=1
       loadct,33,/silent
       for m=0,n_elements(irarray[0,*])-1 do oplot,wave,irarray[*,m],color=m*5
       loadct,0,/silent
       oplot,irout[*,4],irout[*,0]-irout[*,1],linestyle=2
       oplot,irout[*,4],irout[*,0]+irout[*,1],linestyle=2
       oplot,irout[*,4],irout[*,0],thick=2
       oplot,irout[*,4],irout[*,2],thick=2
       oplot,irout[*,4],irout[*,3],thick=2
       wait,0.5
    endif

endfor

stop
END
