PRO plottemplatestars


ans=''
starlist = [0,1,2,3,4] ; This is an array indicating which stars you want plotted

wave1 = 3530.2
disp = 1.125
blue = 3600
red = 5400
;wave1 = 3465.0
;disp = 0.4
;blue = 3300
;red = 6200

files = findfile('vp*cc.fits',count=num)
;stars = fltarr(7000,num)
stars = fltarr(2048,num)
FOR j=0,num-1 DO stars[*,j] = readfits(files[j])

readcol,'Tlist',FORMAT='A,A,A',names,types,metal

wave = fltarr(7000)
FOR j=0,6999 DO wave[j] = wave1 + (j*disp)

window,0,retain=2,xsize=900,ysize=600
device,decomposed=0
loadct,0
REPEAT BEGIN
    w1 = uint((blue-wave1)/disp)
    w2 = uint(((red-wave1)/disp)-1)
    plot,wave[w1:w2],stars[w1:w2,0],xtitle='Wavelength (A)',ytitle='Flux',title='Template Star Spectra',$
      yrange=[0,6],ystyle=1,xrange=[blue-200,red],xstyle=1,/nodata
    FOR j=0,4 DO BEGIN
        indy = starlist[j]
        temp = stars[*,indy]+(j*1.0)
        oplot,wave[w1:w2],temp[w1:w2]
        xyouts,blue-140,j+1,names[indy]
        xyouts,blue-140,j+0.7,'('+types[indy]+')'
    ENDFOR
    print,'Enter w,n or d:'
    read,ans
    IF (ans EQ 'w') THEN BEGIN
        print,'Enter a new blue wavelength:'
        read,blue
        print,'Enter a new red wavelength:'
        read,red
    ENDIF

    IF (ans EQ 'n') THEN BEGIN
        print,'Enter a new set of 5 numbers:'
        read,starlist
    ENDIF
    
ENDREP UNTIL (ans EQ 'd')

print,'The next ENTER deletes the window.'
pause
wdelete,0
stop
end
