pro addspec

;this routine looks for all the #ex.fits files and just adds them
;together.
;it lets you pick the output name

disp1 = 3530
disp2 = 1.125

files = findfile('*ex.fits')
print,transpose(files)

n1 = n_elements(files)
test = readfits(files[0])
n2 = n_elements(test) 
data = dblarr(n2,n1)

wave = fltarr(n2)
for j=0,n2-1 do wave[j] = disp1 + (disp2*j)

for j=0,n1-1 do data[*,j] = readfits(files[j])

ttl = total(data,2)

window,0,retain=2

loadct,0
plot,wave,ttl,xtitle='Wavelength',ytitle='Counts'
loadct,4
for j=0,n1-1 do begin
    oplot,wave,data[*,j],color=60+(j*20)
endfor
loadct,0

print,'Next ENTER deletes the plot.'
pause
wdelete,0

ans=''
print,'Name the output file (w/o the .fits):'
read,ans

name = ans+'.fits'
writefits,name,ttl

end
