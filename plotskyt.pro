; This routine uses the output of "collapse_sky.pro" to overplot the
; subtracted sky values on top of a spectra. Both are in the observed
; frame, so you see where the sky lines fall in your final spectra.

; SKYLIST: The list of the name_skyC.fits frames you got from the
; collapse_sky routine. They need to be interpolated to the same
; wavelength scale. 

; W1 and W2 are just the starting and ending wavelenght ranges.

; BIN: The name of the spectra you want to overplot.

PRO plotskyt, skylist

w1 = 3600
w2 = 5500

wi = 3550
wd = 0.375
swch = 'off'
scale1 = 1.0
scale2 = 1.0

readcol,skylist,format='a',sl
test = readfits(sl[0])

n1= n_elements(sl)
n2 = 6144

skydata = dblarr(n2,n1)
wave = dblarr(n2)

for j=0,n1-1 do skydata[*,j] = readfits(sl[j])

skydata1 = skydata[*,0:22]
skydata2 = skydata[*,23:*]

for j=0,n2-1 do wave[j] = wi+(j*wd)
m = median(skydata,dimension=2,/even)

mm = median(skydata,dimension=1,/even)
mm1 = mm[0:22]
mm2 = mm[23:41]
nm1 = mm1/median(mm1)
nm2 = mm2/median(mm2)
nskydata1 = skydata1
nskydata2 = skydata2

for j=0,22 do nskydata1[*,j] = skydata1[*,j] / nm1[j]
for j=0,18 do nskydata2[*,j] = skydata2[*,j] / nm2[j]

msky1 = median(nskydata1,dim=2)
msky2 = median(nskydata2,dim=2)

msub = msky1 - msky2

nmsky1 = msky1 - median(msky1[4000:5000])
nmsky2 = msky2 - median(msky1[4000:5000])

set_plot,'ps'
device,filename='may-marchSKY.ps',/color
loadct,0
plot,wave,msub,xtitle='Wavelength (A)',ytitle='May Sky - March Sky (Pixels)',xrange=[w1,w2],$
  yrange=[-10,10],xstyle=1,ystyle=1,thick=2,$
  title='Sky Subtraction',xthick=3,ythick=3,charthick=3,charsize=1.3
loadct,4
    oplot,[4102,4102],[-100,100],thick=2,color=150
    oplot,[4341,4341],[-100,100],thick=2,color=150
    oplot,[4861,4861],[-100,100],thick=2,color=150
loadct,0
oplot,wave,msub

;oplot,wave,nmsky1,color=150
;oplot,wave,nmsky2,color=60
device,/close_file
device,filename='monthlySKY.ps',/color
loadct,0
plot,wave,msub,xtitle='Wavelength (A)',ytitle='May Sky - March Sky (Pixels)',xrange=[w1,w2],$
  yrange=[0,100],xstyle=1,ystyle=1,thick=2,/nodata,$
  title='Sky Subtraction',xthick=3,ythick=3,charthick=3,charsize=1.3
loadct,4
    oplot,[4102,4102],[-100,100],thick=2,color=150
    oplot,[4341,4341],[-100,100],thick=2,color=150
    oplot,[4861,4861],[-100,100],thick=2,color=150
oplot,wave,msky1,color=110,thick=2
oplot,wave,msky2,color=150,thick=2
xyouts,5000,40,'May',color=110,charthick=2,charsize=2
xyouts,5000,32,'March',color=150,charthick=2,charsize=2

;oplot,wave,nmsky1,color=150
;oplot,wave,nmsky2,color=60
device,/close_file

stop
end
