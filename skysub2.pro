;THE ONLY DIFFERENCE BETWEEN THIS AND SKYSUB.PRO IS THAT EACH SKY
;COEFFICIENT IS SEPERATE.
;Also, the sky spectra stay for the entire reduction...

PRO skysub2


;This code works from a list, named skysub.list, that sits in the
;'binned directory.  The session is run from the 'reduc' directory

;It's format is a column in the following form

;redshift
;starting wavelength
;dispersion
;sky1
;sky2
;data

;EX:
;0.012843
;3530.1
;1.125
;jm0084.r1.fits
;jm0086.r1.fits
;jm0085.r1.fits

;The output of the code is a file, named vp####ss.fits, with the 'ss'
;for 'sky subtracted'.  This gets put into the 'completed'
;directory. Two plots are also generated.  Both plots have titles
;starting with the science frame, then the radial binning position
;(r1, r2, etc) then the two sky frames used in the background.  Then,
;one has the subscript '.all.ps' while the other is just '.ps'.  The
;'.ps' is just the final spectra.  The '.all.ps' shows the original
;data, the sky that was subtracted, and the resulting subtracted spectrum.

n1=double(2.0)
n2=double(2.0)

ans = ''
cd,'binned'
list = strarr(1,6)
openr,5,'skysub.list'
readf,5,list
free_lun,5

redshift = float(list[0])
wv1 = float(list[1])
disp = float(list[2])
sky1 = readfits(list[3])
sky2 = readfits(list[4])
data = readfits(list[5])
cd,'../'


IF (n_elements(data) EQ 1024) THEN BEGIN
    obswave = fltarr(1024)
    FOR j=0,1023 DO obswave[j]=wv1+(disp*j)
    restwave = obswave-(obswave*redshift)
    dp = 1024
ENDIF
IF (n_elements(data) EQ 2048) THEN BEGIN
    obswave = fltarr(2048)
    FOR j=0,2047 DO obswave[j]=wv1+(disp*j)
    restwave = obswave-(obswave*redshift)
    dp = 2048 
ENDIF

sky = dblarr(dp)
data1 = dblarr(dp)
slope = fltarr(dp)

blue=3600.000
red=5400.000

;The first plot is put up.

mag = 0.0
s = 0.0
sd = s+1.0
sdold=sd
REPEAT BEGIN
        
    FOR j=0,dp-1 DO sky[j] = n1*sky1[j] + n2*sky2[j]
    FOR j=0,dp-1 DO slope[j] = 1 + (j*s)/(dp-1)
    sky = slope*sky
    FOR j=0,dp-1 DO data1[j] = data[j] - sky[j]

    window,3,xsize=900,ysize=450,retain=2
    device,decomposed=0
    loadct,0
    IF (ans NE 'n') THEN BEGIN
        plot,restwave,data1,xrange=[blue,red],title='Sky-subtracted data',$
          xtitle='REST Wavelength',xstyle=1,/ynozero,charsize=1.5,ytitle='CCD Counts'
        xyouts,0.4,0.2,'The current sky coefficients are '+strn(n1)+' and '+strn(n2),/normal,charsize=1.5
    ENDIF
    IF (ans EQ 'n') OR (ans EQ 's') THEN BEGIN
        plot,restwave,data1,xrange=[blue,red],title='Sky-subtracted data',$
          xtitle='REST Wavelength',xstyle=1,/ynozero,charsize=1.5,ytitle='CCD Counts'
        xyouts,0.4,0.2,'The current sky coefficients are '+strn(n1)+' and '+strn(n2),/normal,charsize=1.5
        loadct,4
        oplot,restwave,old,thick=0.5,color=150
        IF (ans EQ 'n') THEN xyouts,0.4,0.15,'The OLD sky coefficients are '+strn(nold1)+' and'+strn(nold2),$
          /normal,charsize=1.5,color=150
        IF (ans EQ 's') THEN xyouts,0.4,0.15,'The OLD sky slope was '+strn(sdold),/normal,charsize=1.5,color=150
    ENDIF

    window,2,xsize=900,ysize=520,retain=2
    device,decomposed=0
    
    loadct,0
    plot,obswave,data,xrange=[blue,red],xstyle=1,title='Sky-subtraction in progress on frame '+list[5], $
      xtitle='OBSERVED Wavelength',charsize=1.5,ytitle='CCD Counts'
    loadct,4
    xyouts,0.15,0.9,'SKY = '+strn(sd)+'('+strn(n1)+' * sky!d1!n + '+strn(n2)+' * sky!d2!n)',/normal,charsize=1.5,color=60

    IF (mag EQ 0) THEN oplot,obswave,sky,color=60

    oplot,obswave,data1,color=110
    xyouts,0.15,0.85,'GREEN is the sky-subtracted DATA.',/normal,color=110,charsize=1.5

    IF (mag NE 0) THEN BEGIN
        oplot,obswave,skymag,color=180
        xyouts,0.15,0.75,'ORANGE: A '+strn(mag)+'X magnification of the current SKY.',/normal,charsize=1.5,color=180
    ENDIF

    oplot,obswave,sky1,color=225
    oplot,obswave,sky2,color=255
;    oplot,obswave,sky
    xyouts,0.15,0.8,'The yellow spectra are the orignial sky frames.',charsize=1.5,color=255,/normal
    
    print,'w: change wavelength,   m: magnify sky,   n: change coefficients'
    print,'s: change slope of the sky,   d: done'
    print,'What next? (w,m,n,s or d)'
    read,ans

    IF (ans NE 'w') AND (ans NE 'm') AND (ans NE 'n') AND (ans NE 's') AND (ans NE 'd') THEN BEGIN
        print,'Not an option!  Try again.'
        read,ans
    ENDIF

    IF (ans EQ 'w') THEN BEGIN
        print,'Enter the starting wavelength:'
        read,blue
        IF (blue EQ 1.0) THEN blue = 3600
        print,'Enter the ending wavelength:'
        read,red
        IF (red EQ 1.0) THEN red = 5400
        print,'The wavelength range is now '+strn(blue)+' to '+strn(red)
    ENDIF

    IF (ans EQ 'm') THEN BEGIN
        print,'Enter a magnification factor for the sky:'
        read,mag
        skymag = sky*mag
    ENDIF
    
    IF (ans EQ 'n') THEN BEGIN
        nold1 = n1
        nold2 = n2
        old = data1
        print,'Enter a new N1:'
;        print,'Enter a new N:'
        read,n1
;        n2=n1
        print,'Enter a new N2:'
        read,n2
    ENDIF

    IF (ans EQ 's') THEN BEGIN
        sdold = sd
        print,'Enter a sky slope:'
        read,sd
        s = sd-1.0
        old = data1
    ENDIF

ENDREP UNTIL (ans EQ 'd')

print,'The fits file of Data1 will be written to the COMPLETED directory.'

cd,'completed'
d1 = strsplit(list[5],'.',/extract)
s1 = strsplit(list[3],'.',/extract)
s2 = strsplit(list[4],'.',/extract)
ftsnm = d1[0]+'ss.'+d1[1]+'.fits'
;ftsnm = d1[0]+'ss.'+d1[2]+'.fits' ;changed this line (1 of 3)
writefits,ftsnm,data1
cd,'../'
 
;print,'Save the plots?'
;read,ans
ans='y'
IF (ans EQ 'y') THEN BEGIN
    red = uint(red)
    blue = uint(blue)
    cd,'completed/plots'
;    wset,2
;    pageinfo = pswindow(/landscape)
    set_plot,'ps'
    nm = d1[0]+'_'+d1[1] ;changed this line (2 of 3)
    device,file=nm+'.all.ps',/color
;    device,_extra=pageinfo
    loadct,0
    plot,obswave,data,xrange=[blue,red],title='Raw Data (black) and estimate of sky (blue)', $
      xtitle='OBSERVED Wavelength (A)',xstyle=1,charsize=1,ytitle='CCD Counts',thick=2,xthick=2,ythick=2,charthick=2
    xyouts,0.15,0.9,'SKY = '+strn(sd)+'('+strn(n1)+' * sky!d1!n + '+strn(n2)+' * sky!d2!n)',/normal,charthick=2
;    xyouts,0.15,0.9,'SKY = '+strn(n1)+'*sky1 + '+strn(n2)+'*sky2',/normal
    xyouts,0.15,0.85,'Data: '+d1[0]+'     Sky frames: '+s1[0]+' & '+s2[0],/normal,charthick=2
    xyouts,0.15,0.8,'Black: original data',/normal,charsize=0.8,charthick=2
    xyouts,0.15,0.76,'Green: sky-subtracted data',/normal,charsize=0.8,charthick=2
    xyouts,0.15,0.72,'Blue: final sky spectrum',/normal,charsize=0.8,charthick=2
    xyouts,0.15,0.68,'Yellow: original sky spectra',/normal,charsize=0.8,charthick=2

    loadct,4
    oplot,obswave,sky,color=60,thick=1.5
    oplot,obswave,data1,color=110,thick=1.5
    oplot,obswave,sky1,color=225,thick=1.5
    oplot,obswave,sky2,color=255,thick=1.5

    IF (mag NE 0) THEN BEGIN
        oplot,obswave,skymag,color=180
        loadct,0
        xyouts,0.15,0.64,'Orange: magnified sky spectrum',/normal,charsize=0.8
    ENDIF
    device,/close_file

;    set_plot,'x'
;    wset,3
;    set_plot,'ps'
;    pageinfo = pswindow(/landscape)
    device,file=nm+'.ps',/color
;    device, _extra=pageinfo
    loadct,0
    plot,restwave,data1,xrange=[blue,red],xstyle=1,title='Sky-subtracted Data', xtitle='REST Wavelength (A)',$
      /ynozero,ytitle='CCD Counts',thick=2,xthick=2,ythick=2,charthick=2
    xyouts,0.15,0.90,'SKY = '+strn(sd)+'('+strn(n1)+' * sky!d1!n + '+strn(n2)+' * sky!d2!n)',/normal,charthick=2
;    xyouts,0.15,0.85,'SKY = '+strn(n1)+'*sky1 + '+strn(n2)+'*sky2',/normal
    xyouts,0.15,0.85,'Sky frames: '+s1[0]+' & '+s2[0],/normal,charthick=2
    xyouts,0.15,0.80,'Data: '+d1[0]+' at '+d1[1],/normal,charthick=2 ;changed this line (3 0f 3)

    device,/close_file
    set_plot,'x'
    cd,'../../'
ENDIF ELSE print,'No plot saved.'

    wdelete,2,3
stop
END

