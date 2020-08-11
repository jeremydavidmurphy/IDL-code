; This routine generates a quick plot of a single 1-D fits
PRO plot1, fits

wi = 3545.0
disp = 1.125

data = readfits(fits)

n1 = n_elements(data)
wave = fltarr(n1)

for j=0,n1-1 do wave[j] = wi + (j * disp)

window,1,retain=2,ypos=50
device,decomposed = 0

x1 = 3700
x2 = 5500
ans = ''
repeat begin
plot,wave,data,xrange=[x1,x2],/xstyle,xtitle='Wavelength (A)',$
  ytitle='Flux',title=fits,/ynozero

print,'Change the plotting range? (y/n)'
read,ans
if (ans eq 'y') then begin
    print,'Enter a new lower and upper wavelength range (follow each by ENTER):'
    read,x1
    read,x2
endif

endrep until (ans eq 'n')

wdelete
stop
END
