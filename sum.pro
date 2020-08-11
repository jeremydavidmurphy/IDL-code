pro sum

;this is a weighted sum of a set of 1-D spectra.
;THEY ARE READ IN FROM A LIST CALLED SUM.LIST

d1 = 3530
d2 = 1.125

readcol,'sum.list',format='A',files

titlefiles = transpose(files)
test = readfits(files[0])
n1 = n_elements(test)
n2 = n_elements(files)

name = ''
ans = ''

wave = fltarr(n1)
data = dblarr(n1,n2)

for j=0,n1-1 do wave[j] = d1 + (d2*j)
for j=0,n2-1 do data[*,j] = readfits(files[j])

wgt = fltarr(n2)
med = fltarr(n2)
out = dblarr(n1)
temp = dblarr(n2)

if (n1 eq 1024) then begin
    for j=0,n2-1 do begin
        piece = data[700:899,j]
        med[j] = median(piece,/even)
    endfor
tl = total(med)
for j=0,n2-1 do wgt[j] = med[j]/tl
endif

if (n1 eq 2048) then begin
    for j=0,n2-1 do begin
        piece = data[1700:1899,j]
        med[j] = median(piece,/even)
    endfor
tl = total(med)
for j=0,n2-1 do wgt[j] = med[j]/tl
endif

for k=0,n1-1 do begin
    onewave = data[k,*]
    for j=0,n2-1 do temp[j] = onewave[j]*wgt[j]
    out[k] = total(temp)
endfor

window,0,retain=2
device,decomposed=0
loadct,0


i1 = where(med eq max(med))
plot,wave,data[*,i1],xtitle='Wavelength',ytitle='Counts',$
  title=files[0]

loadct,4
for j=0,n2-1 do begin
    oplot,wave,data[*,j],color=60+(j*40)
    xyouts,0.2,0.85-(j*0.05),'Spectra #'+strn(j+1),$
      charsize=1.5,color=60+(j*40),/normal
endfor

jj=j

loadct,0
oplot,wave,out
xyouts,0.2,0.85-(jj*0.05),'White is the weighted sum',$
  charsize=1.5,/normal

print,'Like the results? (y or n):'
read,ans

if (ans eq 'y') then begin
    print,'Name the output file (w/o the .fits):'
    read,name
    name = name+'.fits'
    writefits,name,out
    print,''
endif else begin
    print,''
    print,'Modify your sum.list and try again...'
    print,''
endelse

print,'Next ENTER deletes the plot...'
pause
wdelete,0

stop
end
