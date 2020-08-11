; This routine uses the output of "collapse_sky.pro" to overplot the
; subtracted sky values on top of a spectra. Both are in the observed
; frame, so you see where the sky lines fall in your final spectra.

; SKYLIST: The list of the name_skyC.fits frames you got from the
; collapse_sky routine. They need to be interpolated to the same
; wavelength scale. 

; W1 and W2 are just the starting and ending wavelenght ranges.

; BIN: The name of the spectra you want to overplot.

PRO plotsky, skylist, w1, w2, bin, Z=redshift

;EX: 'sky.list', 3800, 4200, 'bin1701.fits'

wi = 3550
wd = 0.375
swch = 'off'
scale1 = 1.0
scale2 = 1.0

readcol,skylist,format='a',sl
test = readfits(sl[0])

data = readfits(bin)

n1= n_elements(sl)
n2 = n_elements(test)

skydata = dblarr(n2,n1)
wave = dblarr(n2)

for j=0,n1-1 do skydata[*,j] = readfits(sl[j])
for j=0,n2-1 do wave[j] = wi+(j*wd)
m = median(skydata,dimension=2,/even)
mm = median(skydata,dimension=1,/even)
nm = mm/median(mm)
nskydata = skydata
for j=0,n1-1 do nskydata[*,j] = skydata[*,j] / nm[j]
msky = median(nskydata,dim=2)

var = dblarr(n2)
nvar = var
for j=0,n2-1 do var[j] = variance(skydata[j,*])
for j=0,n2-1 do nvar[j] = variance(nskydata[j,*])

if (n_elements(redshift) ne 0) then begin
    for j=0,n2-1 do wave[j] = wave[j]-(wave[j]*redshift)
endif
set_plot,'x'
window,0,retain=2,xsize=600,ysize=400,xpos=20,ypos=60
device,decomposed=0
window,1,retain=2,xsize=600,ysize=400,xpos=650,ypos=60
device,decomposed=0

pix1 = round((w1-wi)/wd)
pix2 = round((w2-wi)/wd)
tt = max(skydata[pix1:pix2,*])
ttt = max(data[pix1:pix2])
if (tt gt ttt) then yup=tt else yup=ttt
tt = min(skydata[pix1:pix2,*])
ttt = min(data[pix1:pix2])
if (tt lt ttt) then ydown=tt else ydown=ttt

jump1:
wset,0
loadct,0
datas1 = data * scale1
datas2 = data * scale2

plot,wave,datas1,xtitle='Wavelength (A)',ytitle='Normalized Flux',xrange=[w1,w2],$
  yrange=[ydown,yup],xstyle=1,ystyle=1,thick=1,/nodata,$
  title='WHITE: data  RED: median of sky  GREEN: individual sky frames'

loadct,4
for j=0,n1-1 do begin
    oplot,wave,skydata[*,j],color=110,thick=0.5
endfor
oplot,wave,m,thick=1,color=150
loadct,0
oplot,wave,datas1,thick=1

wset,1
if (swch ne 'on') then begin
    plot,wave,var,xrange=[w1,w2],xstyle=1,title='Variance',xtitle='Wavelength (A)'
    oplot,wave,datas2
    pause
    plot,wave,nvar,xrange=[w1,w2],xstyle=1,title='Normalized Variance',xtitle='Wavelength (A)'
    oplot,wave,datas2
    wset,0
    plot,wave,datas1,xtitle='Wavelength (A)',ytitle='Normalized Flux',xrange=[w1,w2],$
      yrange=[ydown,yup],xstyle=1,ystyle=1,thick=1,/nodata,$
      title='WHITE: data  RED: median of sky  GREEN: individual sky frames'
    loadct,4
    for j=0,n1-1 do begin
        oplot,wave,nskydata[*,j],color=110,thick=0.5
    endfor
    oplot,wave,m,thick=1,color=150
    loadct,0
    oplot,wave,datas1,thick=1
endif

if (swch eq 'on') then begin
    wset,1
    plot,wave,var,xrange=[w1,w2],xstyle=1,yrange=[ydownv,yupv],ystyle=1,$
      title='Variance',xtitle='Wavelength (A)'
    oplot,wave,datas2
    pause
    plot,wave,nvar,xrange=[w1,w2],xstyle=1,title='Normalized Variance',$
      xtitle='Wavelength (A)',yrange=[ydownv,yupv]
    oplot,wave,datas2
    wset,0
    plot,wave,datas1,xtitle='Wavelength (A)',ytitle='Normalized Flux',xrange=[w1,w2],$
      yrange=[ydown,yup],xstyle=1,ystyle=1,thick=1,/nodata,$
      title='WHITE: data  RED: median of sky  GREEN: individual sky frames'
    loadct,4
    for j=0,n1-1 do begin
        oplot,wave,nskydata[*,j],color=110,thick=0.5
    endfor
    oplot,wave,m,thick=1,color=150
    loadct,0
    oplot,wave,datas1,thick=1
endif

ans=''
print,'"w" to change the wavelength range'
print,'"yv" to change the y-variance range.'
print,'"y" to change the y-spectral range.'
print,'"s" to scale the data.'
print,'An ENTER exits the routine.'
read,ans

if (ans eq 'w') then begin
    print,'Wavelength 1:'
    read,w1
    print,'Wavelength 2:'
    read,w2
    goto,jump1
endif
if(ans eq 'yv') then begin
    print,'Yupv:'
    read,yupv
    print,'Ydownv:'
    read,ydownv
    swch = 'on'
    goto,jump1
endif
if (ans eq 'y') then begin
    print,'Yup'
    read,yup
    print,'Ydown'
    read,ydown
    goto,jump1
endif
if (ans eq 's') then begin
    print,'Scale data figure:'
    read,scale1
    print,'Scale variance figure:'
    read,scale2
    goto,jump1
endif

ans=''
print,'Save the plots?'
read,ans
if (ans eq 'y') then begin
    set_plot,'ps'
    device,filename='skyplot.ps',/color
    loadct,0
    plot,wave,data,xtitle='Wavelength (A)',ytitle='Flux (Pixels)',xrange=[w1,w2],$
      yrange=[ydown,yup],xstyle=1,ystyle=1,thick=2,/nodata,$
      title='Sky Subtraction',xthick=3,ythick=3,charthick=3,charsize=1.3
    for j=0,n1-1 do begin
        oplot,wave,skydata[*,j],color=200,thick=0.5
    endfor
    loadct,4
    oplot,wave,msky,thick=3,color=150
    loadct,0
    oplot,wave,datas1,thick=2
;    oplot,[4102,4102],[0,100],thick=2
;    oplot,[4341,4341],[0,100],thick=2
;    oplot,[4861,4861],[0,100],thick=2
    device,/close_file
    
    device,filename='sky_variance.ps',/color
    loadct,0
        if (swch eq 'off') then begin
            loadct,0,/silent
            plot,wave,nvar,xrange=[w1,w2],xstyle=1,title='Normalized Sky Variance',xtitle='Observed Wavelength (A)',$
          charthick=3,xthick=3,ythick=3,thick=2,charsize=1.3
            loadct,4,/silent
            oplot,wave,datas2,color=150
    endif
    if (swch eq 'on') then begin
        loadct,0,/silent
        plot,wave,nvar,xrange=[w1,w2],xstyle=1,yrange=[ydownv,yupv],ystyle=1,$
          title='Normalized Sky Variance',xtitle='Observed Wavelength (A)',charthick=3,xthick=3,ythick=3,thick=2
        loadct,4,/silent
        oplot,wave,datas2,color=150
    endif
    device,/close_file
endif
set_plot,'x'
pause
wdelete,0,1

free_lun,5
openw,5,'M87f_variance.txt'
for j=0,n_elements(var)-1 do printf,5,wave[j],nvar[j]
free_lun,5
stop
end
