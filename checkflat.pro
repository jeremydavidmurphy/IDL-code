pro checkflat,x1, x2, nw

bias=mrdfits('bias.fits', 0, h)

flat0=mrdfits('test0001.fits', 0, h); used to nomalize all the others

flat1=mrdfits('test0022.fits', 0, h) ; blue
flat2=mrdfits('test0026.fits', 0, h) ; cyan
flat3=mrdfits('test0035.fits', 0, h) ; green
flat4=mrdfits('test0038.fits', 0, h) ; orange
flat5=mrdfits('test0040.fits', 0, h) ; red


flat0=flat0-bias
flat1=flat1-bias
flat2=flat2-bias
flat3=flat3-bias
flat4=flat4-bias
flat5=flat5-bias

yprof0=fltarr(2048)
yprof1=fltarr(2048)
yprof2=fltarr(2048)
yprof3=fltarr(2048)
yprof4=fltarr(2048)
yprof5=fltarr(2048)

;x1=500
;x2=520

for i=0, 2048-1 do yprof0[i]=median(flat0[x1:x2,i])
for i=0, 2048-1 do yprof1[i]=median(flat1[x1:x2,i])
for i=0, 2048-1 do yprof2[i]=median(flat2[x1:x2,i])
for i=0, 2048-1 do yprof3[i]=median(flat3[x1:x2,i])
for i=0, 2048-1 do yprof4[i]=median(flat4[x1:x2,i])
for i=0, 2048-1 do yprof5[i]=median(flat5[x1:x2,i])


y=findgen(2048)

aux0=yprof0[sort(yprof0)]
aux1=yprof1[sort(yprof1)]
aux2=yprof2[sort(yprof2)]
aux3=yprof3[sort(yprof3)]
aux4=yprof4[sort(yprof4)]
aux5=yprof5[sort(yprof5)]

norm0=2.0*median(aux0[0:2047])
norm1=2.0*median(aux1[0:2047])
norm2=2.0*median(aux2[0:2047])
norm3=2.0*median(aux3[0:2047])
norm4=2.0*median(aux4[0:2047])
norm5=2.0*median(aux5[0:2047])

;norm0=total(flat0)/total(flat0)*median(yprof0)*2.0
;norm1=total(flat1)/total(flat0)*median(yprof0)*2.0
;norm2=total(flat2)/total(flat0)*median(yprof0)*2.0
;norm3=total(flat3)/total(flat0)*median(yprof0)*2.0
;norm4=total(flat4)/total(flat0)*median(yprof0)*2.0
;norm5=total(flat5)/total(flat0)*median(yprof0)*2.0

;window, nw, retain=2, xsize=1200, ysize=600
device, decompose=0
loadct, 39

set_plot,'ps'
device,file='profile2.ps', xsize=24, ysize=15, xoffset=0, yoffset=0, /color

plot, y, yprof1/norm1, xrange=[0,700], position=[0.03,0.73,0.98,0.99], /norm, xstyle=1, yrange=[-0.1,1.1], ystyle=1
oplot, y, yprof1/norm1, color=50
oplot, y, yprof2/norm2, color=100
oplot, y, yprof3/norm3, color=150
oplot, y, yprof4/norm4, color=200
oplot, y, yprof5/norm5, color=250

plot, y, yprof1/norm1, xrange=[700,1400], position=[0.03,0.38,0.98,0.64], /norm, /noerase, xstyle=1, yrange=[-0.1,1.1], ystyle=1
oplot, y, yprof1/norm1, color=50
oplot, y, yprof2/norm2, color=100
oplot, y, yprof3/norm3, color=150
oplot, y, yprof4/norm4, color=200
oplot, y, yprof5/norm5, color=250

plot, y, yprof1/norm1, xrange=[1400,2048], position=[0.03,0.04,0.98,0.29], /norm, /noerase, xstyle=1, yrange=[-0.1,1.1], ystyle=1
oplot, y, yprof1/norm1, color=50
oplot, y, yprof2/norm2, color=100
oplot, y, yprof3/norm3, color=150
oplot, y, yprof4/norm4, color=200
oplot, y, yprof5/norm5, color=250

device, /close_file
set_plot, 'x'

set_plot,'ps'
device,file='profile2_zooms.ps', xsize=24, ysize=15, xoffset=0, yoffset=0, /color

plot, y, yprof1/norm1, xrange=[0,400], position=[0.03,0.73,0.98,0.99], /norm, xstyle=1, yrange=[-0.1,1.1], ystyle=1
oplot, y, yprof1/norm1, color=50
oplot, y, yprof2/norm2, color=100
oplot, y, yprof3/norm3, color=150
oplot, y, yprof4/norm4, color=200
oplot, y, yprof5/norm5, color=250

plot, y, yprof1/norm1, xrange=[900,1300], position=[0.03,0.38,0.98,0.64], /norm, /noerase, xstyle=1, yrange=[-0.1,1.1], ystyle=1
oplot, y, yprof1/norm1, color=50
oplot, y, yprof2/norm2, color=100
oplot, y, yprof3/norm3, color=150
oplot, y, yprof4/norm4, color=200
oplot, y, yprof5/norm5, color=250

plot, y, yprof1/norm1, xrange=[1648,2048], position=[0.03,0.04,0.98,0.29], /norm, /noerase, xstyle=1, yrange=[-0.1,1.1], ystyle=1
oplot, y, yprof1/norm1, color=50
oplot, y, yprof2/norm2, color=100
oplot, y, yprof3/norm3, color=150
oplot, y, yprof4/norm4, color=200
oplot, y, yprof5/norm5, color=250

device, /close_file
set_plot, 'x'


stop
end

