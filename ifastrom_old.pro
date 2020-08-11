pro ifastrom, list, yindex
; this routine reduces the imaging fiber astrometry data.
;the list is the list of frames. The code will look for a file named
;'coords.txt' that lists the input X and Y position. This is used to
;then associate the input R to the output R.
;THE COORDS.TXT FILE MUST BE IN THE SAME ORDER AS THE FILES IN LIST.

;The routine will look for 'yindex.txt', which is just the input y
;positions, each one listed just once.

;These values set the region on the chip over which to look for a spot.
;test 1 range
;rx1 = 1220
;rx2 = 1600
;ry1 = 1070
;ry2 = 1500
;test 2 range
rx1 = 1224
rx2 = 1579
ry1 = 1074
ry2 = 1449
;test 3 range
;rx1 = 1360
;rx2 = 1540
;ry1 = 1200
;ry2 = 1350

readcol,'coords.txt',format='i,i,i,i',Xind,Yind,Xin,Yin
;readcol,'yindex.txt',format='i',yindex
readcol,list,format='a',files

n1 = n_elements(files)
;n2 = n_elements(yindex)

specs = fltarr(7,n1)
dataarr = fltarr(abs(rx1-rx2)+1,abs(ry1-ry2)+1,n1)

;the gaussian centroids are found
for j=0,n1-1 do begin
    temp = readfits(files[j],/silent)
    trim = temp[rx1:rx2,ry1:ry2]
    dataarr[*,*,j] = trim
    print,'Fitting file '+files[j]
    out = gauss2dfit(trim,ts,/tilt)
    specs[*,j] = ts
endfor

;pixels are converted into microns
specs[2,*] = specs[2,*] * 12.0
specs[3,*] = specs[3,*] * 12.0
specs[4,*] = specs[4,*] * 12.0
specs[5,*] = specs[5,*] * 12.0

;a rough composite frame is created, after a rough background is
;subtracted
dataarr = dataarr - min(dataarr)
dataout = max(dataarr,dimension=3)
writefits,list+'.fits',dataout

;a first attempt at rotation and magnification is made...
;NOTE: THIS IS DEPENDENT ON THE ORDER THE DATA IS TAKEN AND CURRENT
;(FEB 1, 2009) ORIENTATION OF THE IMAGING BUNDLE!!!! THE X AND Y INDEX
;THAT BECOMES THE POINT OF ROTATION MAY CHANGE!

xout = specs[4,*]
yout = specs[5,*]

thetarr = 0

;an estimate of the degree of rotation is made
for j=0,n2-1 do begin
    rowi = where(Yin eq yindex[j])
    n3 = n_elements(rowi)
    if (n3 ge 2) then begin
        for k=0,n3-1 do begin
            sxout = xout[rowi] - xout[rowi[k]]
            syout = yout[rowi] - yout[rowi[k]]
            nozeroi = where(sxout ne 0.0)
            for l=0,n_elements(nozeroi)-1 do begin
                ind = nozeroi[l]
                thetarr = [thetarr,atan(abs(sxout[ind]/syout[ind]))]
            endfor
        endfor
    endif
endfor

thetarr = thetarr[1:*]
theta = median(thetarr)*57.2957795
print,theta
pause

;the rotation is now set.

;the 'results' array is a dx, dy and resulting dR array calculated by
;subtracting every x and y position in the 'specs' array from itself

radius = fltarr(n1-1)
for j=0,n1-1 do begin
   radius[j] = sqrt(xout[j]^2+yout[j]^2)
endfor

stop

results = fltarr(6,n1,n1)

for j=0,n1-1 do begin
    for k=0,n1-1 do begin
        results[0,j,k] = abs(Xin[j] - Xin[k]);dX_in
        results[1,j,k] = abs(Yin[j] - Yin[k]);dY_in
        results[2,j,k] = sqrt(results[0,j,k]^2 + results[1,j,k]^2);dR_in
        results[3,j,k] = abs(specs[4,j] - specs[4,k]);dX_out
        results[4,j,k] = abs(specs[5,j] - specs[5,k]);dY_out
        results[5,j,k] = sqrt(results[3,j,k]^2 + results[4,j,k]^2);dR_out
    endfor
endfor




;This x-index (in this case, the second zero) sets the spot image
;about which you are going to rotate your frame and align with the 0,0
;position of the input fiber. This is what will change if your
;zero-point changes.
rspecs[0,*] = rspecs[0,*] - rspecs[0,0]
rspecs[1,*] = rspecs[1,*] - rspecs[1,0]

;an estimate of the required rotation is calculated, based on
;individual measures of the angles for the first row of spots
theta = fltarr(4)
for j=1,4 do begin
    theta[j-1] = atan(rspecs[1,j]/rspecs[0,j]) * 57.295779
endfor
print,theta
theta = mean(theta)
rtheta = theta * (-1.0)
;the magnification of the system is estimated, again by the average
;offset between the spots on the first row.

mag = fltarr(4)
for j=1,4 do begin
    mag[j-1] = rspecs[0,j]/Xin[j]
endfor
print,mag
mag = mean(mag)
rmag = 1/mag

rotarr = dataarr
ix = where(specs[0,*] eq specs[0,0])
for j=0,n1-1 do rotarr[*,*,j] = rot(dataarr[*,*,j],rtheta,rmag,specs[0,0],$
                                    specs[1,0],/pivot,/interp)
dataoutr = max(rotarr,dimension=3)
writefits,list+'rot.fits',dataoutr
stop


;the radius array is sorted. the zeros are tossed
Rin = results[2,*,*]
Rout = results[5,*,*]
Rins = Rin[bsort(Rin)]
Routs = Rout[bsort(Rin)]
index = where(Routs ne 0.0)
Rins = Rins[index]
Routs = Routs[index]

n4 = n_elements(Rins)

;the repeats in the array (due to its symmetry) are tossed. They are
;also placed into a single array
Rall = fltarr(2,n4/2.0)
n5 = n_elements(Rall[0,*])
cntr = 0
for j=0,n4-1,2 do begin
    Rall[0,cntr] = Rins[j]
    Rall[1,cntr] = Routs[j]
    cntr = cntr + 1
endfor

;The input and output radius measurements are written to a text file.
openw,5,list+'.rad'
for j=0,n4-1 do printf,5,Rins[j],Routs[j]
free_lun,5

;the size of the output stats is determined
cntr = 0
for j=0,n5-2 do begin
    ele1 = Rall[0,j]
    ele2 = Rall[0,j+1]
    if (ele1 ne ele2) then cntr = cntr + 1
endfor
outarr0 = fltarr(6,cntr)
outarr1 = fltarr(6,cntr)
n4 = cntr
;the index for the steps is calculated
cntr1 = 0
cntr2 = 0
cntr3 = 0
while (cntr1 lt n5) do begin
    repeat begin
        if (cntr1+cntr2+1 lt n5) then begin
            delta = Rall[0,cntr1+cntr2+1] / Rall[0,cntr1+cntr2]
            cntr2 = cntr2 + 1
        endif else begin
            if cntr2-1 lt 2 then goto,jump1
        endelse
    endrep until (delta ne 1.00)
    piece0 = Rall[0,cntr1:cntr1+cntr2-1];the input R-values
    piece1 = Rall[1,cntr1:cntr1+cntr2-1];the resulting R-values
    outarr0[*,cntr3] = [n_elements(piece0),mean(piece0),median(piece0),$
                       stddev(piece0),min(piece0),max(piece0)]
    outarr1[*,cntr3] = [n_elements(piece1),mean(piece1),median(piece1),$
                       stddev(piece1),min(piece1),max(piece1)]
    cntr1 = cntr1+cntr2
    cntr2 = 0
    cntr3 = cntr3 + 1
    
endwhile
jump1:

index = bsort(outarr0[1,*])
outarr0s = outarr0[*,index]
outarr1s = outarr1[*,index]
f0 = '(a11,a10,a10,a13,a10,a12,a12)'
f1 = '(i3,8x,f10.4,f10.4,3x,f10.5,2x,f8.5,2x,f10.5,2x,f10.5)'
openw,5,list+'.stats'
printf,5,['N_Elements','Input Rad','Mean Rad','Median Rad','Std Dev','Min','Max'],format=f0
for j=0,n4-1 do begin
;    printf,5,outarr0s[*,j],format=f1
;    printf,5,outarr1s[*,j],format=f1
    printf,5,outarr1s[0,j],outarr0s[1,j],outarr1s[1:*,j],format=f1
endfor
free_lun,5

coeffs = linfit(outarr0s[1,*],outarr1s[1,*],chisq=chi,yfit=yvals)
;plots are generated...
window,0,retain=2
device,decomposed=0

plot,rall[0,*],rall[1,*],psym=3,/ynozero
oplot,outarr0s[1,*],yvals
;xyouts,outarr0s[1,0],(max(outarr0s[1,*])+0.5*max(outarr0s[1,*])),'coeffs','chi-squared'
xyouts,outarr0s[1,0],(max(outarr0s[1,*])-0.15*max(outarr0s[1,*])),coeffs[0]
xyouts,outarr0s[1,0],(max(outarr0s[1,*])-0.18*max(outarr0s[1,*])),coeffs[1]
xyouts,outarr0s[1,0],(max(outarr0s[1,*])-0.21*max(outarr0s[1,*])),chi

;a plotting attempt to undo the magnification and rotation of the
;output.

window,0,retain=2
device,decomposed=0

angle = 0.0
mag = 1.0

ans = ''
ans2 = ''

outpoints = fltarr(2,n_elements(specs[4,*]))
outpoints[0,*] = specs[4,*] - min(specs[4,*])
outpoints[1,*] = specs[5,*] - min(specs[5,*])
xoff = 0.0
yoff = 0.0
repeat begin
loadct,0
plot,Xin,Yin,psym=2,xrange=[min(xin)-50,max(xin)+50],xstyle=1,$
  yrange=[min(yin)-50,max(yin)+50],ystyle=1

rotspec = rot(outpoints,angle,mag)
loadct,4
oplot,rotspec[0,*]+xoff,rotspec[1,*]+yoff,psym=4,color=150
print,''
print,'The rotated X and Y min values are:'
print,min(rotspec[0,*])
print,min(rotspec[1,*])
print,''
print,'Enter [a]ngle, [m]agnification, [o]ffset, [d]one'
read,ans
if (ans eq 'a') then begin
    print,'Enter a new angle.'
    read,angle
endif
if (ans eq 'm') then begin
    print,'Enter a new magnification.'
    read,mag
endif
if (ans eq 'o') then begin
    print,'[x] or [y] or [b]oth?'
    read,ans2
    if (ans2 eq 'x') then begin
        print,'Enter a new x-offset'
        read,xoff
    endif
    if (ans2 eq 'y') then begin
        print,'Enter a new y-offset'
        read,yoff
    endif
    if (ans2 eq 'b') then begin
        print,'Enter a new x-offset'
        read,xoff
        print,'Enter a new y-offset'
        read,yoff
    endif
endif

endrep until (ans eq 'd')

stop
end
