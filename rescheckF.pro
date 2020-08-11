; This function reads in a name of the form:
; _jan08_n2
; It then reads in the corresponding arc, ptow and mask file and
; calculates the FWHM of the lines in the HgCd.dat file (these are now
; internal so as to avoid confusion) and returns both FWHM and
; dispersion in km/sec as a n,m,2 array
; n: wavelength (11 of them- see below for their values)
; m: fiber number, with fiber #1 being at the top
; 2: FWHM and dispersion

FUNCTION rescheckF, name

splot = 'on' ;switch this to 'on' to plot each fit to screen. 
splot = 'off' ;switch this to 'on' to plot each fit to screen. 

nprofile = 5 ;the number of rows in a fiber profile
step = (nprofile-1)/2

hgcd1 = [4046.5539, 4077.8298, 4358.3253, 4678.1474,$
         4799.9080, 5085.8173, 5460.7366, 5769.6000, 5790.6239]
hgcd1 = [4046.5539, 4077.8298, 4358.3253,4799.9080, 5085.8173, 5769.6]

temp = strsplit(name,'.',/extract)
name = temp[0]

arc = 'arc_'+name+'e.fits'
ptow = 'ptow_'+name+'.dat'
mask = 'mask_'+name+'.dat'

data = readfits(arc,/silent)
nf0 = n_elements(data[*,0])

readcol,mask,silent=1,format='i,i',row,fib
nf1 = n_elements(fib)

index = row - 1

readcol,ptow,silent=1,format='f,f,f,f,f',c0,c1,c2,c3,c4
if (n_elements(c4) ne nf1) then stop

wave = fltarr(nf0,nf1)
for jj=0,nf1-1 do begin
    for k=0,nf0-1 do begin
        wave[k,jj] = c0[jj]+c1[jj]*k+c2[jj]*k*k+c3[jj]*k*k*k+c4[jj]*k*k*k*k
    endfor
endfor
disp = median(c1)

nf2 = n_elements(hgcd1)
hgcd0 = hgcd1-12
hgcd2 = hgcd1+12

collapsed = fltarr(nf0,nf1)
values = fltarr(nf2,nf1,4) ; wavelength,fiber #,(wave,deltawave,FWHM,disp)

print,'Fitting gaussians to '+arc
if (splot eq 'on') then window,2,retain=2

for jj=0,nf1-1 do begin ;a loop through fibers
    fiber = total(data[*,index[jj]-step:index[jj]+step],2)
    onew = wave[*,jj]
    if (total(onew) eq 0) then goto, jump1
    for kk=0,nf2-1 do begin ;a loop through each wavelength in the arcs
        i = where(onew gt hgcd0[kk] and onew lt hgcd2[kk])
        pd = fiber[i]
        pw = onew[i]
        temp = gaussfit(pw,pd,A,nterms=4)
        if (splot eq 'on') then begin
            plot,pw,pd,thick=2
            oplot,pw,temp,linestyle=3,thick=2
;            wait,0.05
        endif
        fwhmA = 2*sqrt(2*alog(2))*A[2]
        values[kk,jj,0] = A[1]
        values[kk,jj,1] = A[1]-hgcd1[kk]
        values[kk,jj,2] = fwhmA
        values[kk,jj,3] = 299792.458*(A[2]*disp/A[1])
    endfor
    jump1: ;skips over dead fibers
endfor

i = where(values[0,*,0] eq 0)
if (i[0] ne -1) then values[*,i,*] = !Values.F_NAN ;the dead fibers are set to NaN's

;out = values[*,*,2:3] ;just the FWHM and dispersion are returned
out = values
return,out

END
