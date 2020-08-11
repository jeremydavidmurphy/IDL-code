;This routine was used to determine the astrometry quality of the
;imaging fiber bundle.

pro ifastrom, list
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
;radius = fltarr(5,5,25)
;nnn=4

;test 2 range
rx1 = 1224
rx2 = 1579
ry1 = 1074
ry2 = 1449
radius = fltarr(5,5,25)
nnn=4

;test 3 range
;rx1 = 1360
;rx2 = 1540
;ry1 = 1200
;ry2 = 1350
;radius = fltarr(10,10,100)
;nnn=9

;readcol,'coords.txt',format='i,i,i,i',xi,yi,xin,yin
readcol,list,format='a',files
n1 = n_elements(files)

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

stop
xout = specs[4,*]
yout = specs[5,*]

;test1
xgrid = [[-1,xout[0:1],-1,-1],[xout[2:5],-1],[xout[6:10]],[xout[11:14],-1],[-1,xout[15],-1,-1,-1]]
ygrid = [[-1,yout[0:1],-1,-1],[yout[2:5],-1],[yout[6:10]],[yout[11:14],-1],[-1,yout[15],-1,-1,-1]]


;test2
;xgrid = [[xout[0:4]],[xout[5:9]],[xout[10:14]],[xout[15:19]],[xout[20:24]]]
;ygrid = [[yout[0:4]],[yout[5:9]],[yout[10:14]],[yout[15:19]],[yout[20:24]]]

;test3
;xgrid =[[xout[0:9]],[xout[10:19]],[xout[20:29]],[xout[30:39]],[xout[40:49]],$
;        [xout[50:59]],[xout[60:69]],[xout[70:79]],[xout[80:89]],[xout[90:99]]]
;ygrid =[[yout[0:9]],[yout[10:19]],[yout[20:29]],[yout[30:39]],[yout[40:49]],$
;        [yout[50:59]],[yout[60:69]],[yout[70:79]],[yout[80:89]],[yout[90:99]]]

;a rough composite frame is created, after a rough background is
;subtracted
dataarr = dataarr - min(dataarr)
dataout = max(dataarr,dimension=3)
writefits,list+'.fits',dataout

l=0
for k=0,nnn do begin
    for j=0,nnn do begin
        xo = xgrid[j,k]
        yo = ygrid[j,k]
        if (xo eq -1.0) then goto, jump1
        for i=0,nnn do begin
            xi = xgrid[i,k]
            yi = ygrid[i,k]
            if (xi eq -1.0) then goto, jump2
            radius[i,k,l] = sqrt((xo-xi)^2 + (yo-yi)^2)
            jump2:
            xi = xgrid[j,i]
            yi = ygrid[j,i]
            if (xi eq -1.0) then goto, jump3
            radius[j,i,l] = sqrt((xo-xi)^2 + (yo-yi)^2)
            jump3:
        endfor
        jump1:
        l = l + 1
    endfor
endfor

radius = reverse(radius,2)

;plotting routines for each test...
;test3

;set_plot,'ps'
;device,file='test3.ps',/color

;loadct,0
;test = histogram(radius,locations=loc,min=2,max=95)
;plot,loc,test,psym=10,xrange=[0,100],xstyle=1,xtitle='!6Spot Separation (!7l!6m)',$
;  ytitle='Number',title='Histogram for 10!7l!6m step tests',xthick=2,ythick=2,$
;  charsize=1.2,charthick=2,thick=2

;device,/close_file
;set_plot,'x'


;test1
;this is done with plotter.pro found in the test1 directory

stop
end
