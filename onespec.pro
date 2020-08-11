PRO onespec, F=file, w1=w1, w2=w2, redshift=redshift

;this program plots 1-d fits files.  it can handle up to 5.
;It requires all spectra to be at the same Z and be of the same length.

IF (n_elements(w1) EQ 0) THEN w1 = 3530
IF (n_elements(w2) EQ 0) THEN w2 = 1.125



blue = 3600
red = 5400

cntr=0
ans=''


IF (N_Elements(file) EQ 0) THEN BEGIN
    file = strarr(1)
    print,'Enter the file:'
    read,file
    file = file+'.fits'
ENDIF

IF (N_Elements(redshift) EQ 0) THEN BEGIN
    redshift = fltarr(1)
    print,'Enter the redshift:'
    read,redshift
ENDIF

color=[60,110,150,180,255]

d1 = readfits(file)
n = n_elements(d1)
data = fltarr(n,5)
data[*,0] = d1
names = strarr(5)
names[0] = file
wave = fltarr(n)
restwave = wave
FOR j=0,n-1 DO wave[j] = w1+(w2*j)
FOR j=0,n-1 DO restwave[j] = wave[j]-(wave[j]*redshift)

window,2,retain=2,xsize=900,ysize=600
device,decomposed=0
loadct,0
plot,restwave,d1,xrange=[blue,red],xstyle=1,title='Spectra for file '+file,$
  xtitle='Rest Wavelength (A)',/nodata,/ynozero
loadct,4
oplot,restwave,d1,color=color[0]


REPEAT BEGIN
    print,'Select an option (w, n, y, d):'
    print,'(W:wavelength, N:new spectra, Y:rescale Y-axis, D:Done)'
    read,ans
    IF (ans EQ 'n') THEN BEGIN
        newfile=''
        cntr = cntr+1
        print,'Enter the new file name:'
        read,newfile
        newfile = newfile+'.fits'
        names[cntr] = newfile
        data[*,cntr] = readfits(newfile)
        loadct,0
        IF (n_elements(yup) EQ 0) THEN BEGIN
            plot,restwave,data[*,0],xrange=[blue,red],xstyle=1,$
              xtitle='Rest Wavelength (A)',title='Spectra for file '+file,/nodata,/ynozero
            loadct,4
            FOR k = 0,cntr DO oplot,restwave,data[*,k],color=color[k]
        ENDIF ELSE BEGIN
            plot,restwave,data[*,0],xrange=[blue,red],yrange=[ydown,yup],ystyle=1,xstyle=1,$
              xtitle='Rest Wavelength (A)',title='Spectra for file '+file,/nodata
            loadct,4
            FOR k = 0,cntr DO oplot,restwave,data[*,k],color=color[k]
        ENDELSE
    ENDIF

    IF (ans EQ 'w') THEN BEGIN
        print,'Enter a new blue wavelength:'
        read,blue
        print,'Enter a new red wavelength:'
        read,red
        loadct,0
       IF (n_elements(yup) EQ 0) THEN BEGIN
            plot,restwave,data[*,0],xrange=[blue,red],xstyle=1,$
              xtitle='Rest Wavelength (A)',title='Spectra for file '+file,/nodata,/ynozero
            loadct,4
            FOR k = 0,cntr DO oplot,restwave,data[*,k],color=color[k]
        ENDIF ELSE BEGIN
            plot,restwave,data[*,0],xrange=[blue,red],yrange=[ydown,yup],ystyle=1,xstyle=1,$
              xtitle='Rest Wavelength (A)',title='Spectra for file '+file,/nodata
            loadct,4
            FOR k = 0,cntr DO oplot,restwave,data[*,k],color=color[k]
        ENDELSE
     ENDIF

    IF (ans EQ 'y') THEN BEGIN
        yup = dblarr(1) & ydown = yup
        print,'Enter new upper limit:'
        read,yup
        print,'Enter a new lower limit:'
        read,ydown
        loadct,0
       IF (n_elements(yup) EQ 0) THEN BEGIN
            plot,restwave,data[*,0],xrange=[blue,red],xstyle=1,$
              xtitle='Rest Wavelength (A)',title='Spectra for file '+file,/nodata,/ynozero
            loadct,4
            FOR k = 0,cntr DO oplot,restwave,data[*,k],color=color[k]
        ENDIF ELSE BEGIN
            plot,restwave,data[*,0],xrange=[blue,red],yrange=[ydown,yup],ystyle=1,xstyle=1,$
              xtitle='Rest Wavelength (A)',title='Spectra for file '+file,/nodata
            loadct,4
            FOR k = 0,cntr DO oplot,restwave,data[*,k],color=color[k]
        ENDELSE
    ENDIF

ENDREP UNTIL (ans EQ 'd')

print,'Save the plot?'
read,ans
IF (ans EQ 'y') THEN BEGIN
    ans1=''
    print,'Add a note to the plot?'
    read,ans1
    IF (ans1 EQ 'y') THEN BEGIN
        note = strarr(1)
        print,'Type away...'
        read,note
        note = 'NOTE: '+note
    ENDIF
    IF (cntr EQ 0) THEN t = 'Spectra for file '+file
    IF (cntr NE 0) THEN t = 'Spectra for various files'
    set_plot,'ps'
    device,file='spec.ps',/color
    loadct,0
    IF (n_elements(yup) EQ 0) THEN BEGIN
        plot,restwave,d1,xrange=[blue,red],xstyle=1,title=t,$
          xtitle='Rest Wavelength (A)',/nodata,/ynozero,xthick=2,ythick=2,thick=2,charthick=2
        loadct,4
        FOR k = 0,cntr-1 DO oplot,restwave,data[*,k],color=color[k],thick=2
        FOR k = 0,cntr-1 DO xyouts,0.2,0.90-(0.05*k),names[k],/normal,color=color[k],charthick=2
        IF (cntr EQ 4) THEN BEGIN
            loadct,0
            oplot,restwave,data[*,4],thick=2
            xyouts,0.2,0.7,names[4],/normal,charthick=2
        ENDIF
    ENDIF ELSE BEGIN
        plot,restwave,d1,xrange=[blue,red],yrange=[ydown,yup],ystyle=1,xstyle=1,title=t,$
          xtitle='Rest Wavelength (A)',/nodata,xthick=2,ythick=2,thick=2,charthick=2
        loadct,4
        FOR k = 0,cntr-1 DO oplot,restwave,data[*,k],color=color[k],thick=2
        FOR k = 0,cntr-1 DO xyouts,0.2,0.90-(0.05*k),names[k],/normal,color=color[k],charthick=2
        IF (cntr EQ 4) THEN BEGIN
            loadct,0
            oplot,restwave,data[*,4],thick=2
            xyouts,0.2,0.7,names[4],/normal,charthick=2
        ENDIF
    ENDELSE
    IF (ans1 EQ 'y') THEN BEGIN
        loadct,0
        xyouts,0.15,0.15,note,charsize=0.7,/normal
    ENDIF
    device,/close
    set_plot,'x'
    print,'The plot will be written as spec.ps'
ENDIF ELSE print,'No plot for you!'
wdelete,2
stop

END



