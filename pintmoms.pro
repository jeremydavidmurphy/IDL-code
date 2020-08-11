PRO pintmoms, file, scale

; This code is used to plot the intmom.out file. This is the internal
; moments of the modeling results. The difference between this and
; pintmom.pro (i.e. the 's') is that it plots a single intmom file.

; FILE: The name of the intmom.out file
; SCALE: This is the size scale of the galaxy. EX: M87 has a scale of
; 86.5 pc/arcsec. This allows a second x-axis plotted showing kpc.

x1 = 0.01 ;the minimim plotting range
x2 = 800 ;the maximum plotting range
y1 = 0.2
y2 = 1.8
;********************************************************************

angles = [5.77,17.56,30.22,45.00,71.57]

readcol,file,format='f,f,f,f,x,f,f,f',R,theta,sigR,St,Sp,vp,beta,skipline=1

n0 = n_elements(R)
n1 = n_elements(angles)
n2 = n0 / n1

sigarr = dblarr(n2,2,n1)
radius = fltarr(n2,2)

cntr1 = 0
cntr2 = 0
for j=0,n0-1 do begin
    t1 = theta[j]
    t2 = angles[cntr1]
;    print,t1,t2
    if (t1 eq t2) then begin
        sigarr[cntr2,0,cntr1] = sigR[j]
        sigarr[cntr2,1,cntr1] = sqrt((St[j]^2 + Sp[j]^2 + vp[j]^2)*0.5)
;        print,sigarr[cntr2,0,cntr1]
;        print,sigarr[cntr2,1,cntr1]
        cntr2 = cntr2 + 1
        if (cntr1 eq 0) then radius[j,0] = R[j]
    endif else begin
        cntr2 = 0
        cntr1 = cntr1 + 1
        sigarr[cntr2,0,cntr1] = sigR[j]
        sigarr[cntr2,1,cntr1] = sqrt((St[j]^2 + Sp[j]^2 + vp[j]^2)*0.5)
;        print,sigarr[cntr2,0,cntr1]
;        print,sigarr[cntr2,1,cntr1]
        cntr2 = 1
    endelse
endfor

for j=0,n2-1 do begin
    RR = sigarr[j,0,*]
    TH = sigarr[j,1,*]
    div = RR / TH
    radius[j,1] = mean(div)
endfor

radkpc1 = (x1 * scale) / 1000.0
radkpc2 = (x2 * scale) / 1000.0

set_plot,'ps'
device,file='intmoms.ps'
loadct,0,/silent
plot,radius[*,0],radius[*,1],/xlog,xrange=[x1,x2],xstyle=9,$
  xtitle='!6Radius (arcsec)',ytitle='!7r!3!dr!n/!7r!3!dt',$
  charsize=1.2,xthick=3,ythick=3,thick=5,charthick=3,/nodata,$
  yrange=[y1,y2],/ys,position=[0.13,0.13,0.93,0.89]
loadct,4,/silent
oplot,radius[*,0],radius[*,1],thick=5,color=150
loadct,0,/silent
axis,xaxis=1,xrange=[radkpc1,radkpc2],/save,xtitle='!6R (kpc)',$
  charthick=3,charsize=1.2,xstyle=1,xthick=3
device,/close_file

STOP
END
