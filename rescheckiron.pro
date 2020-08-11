; This version is specific for the calcRES.pro routine.

FUNCTION rescheckiron, arc,ptow,mask,Splot=splot

wavelist = [3592.32,3731.13,4556.23,4620.78,4669.06,4776.03]

if (n_elements(Splot) eq 0) then splot = 'off' else splot='on'

nprofile = 5 ;the number of rows in a fiber profile
step = (nprofile-1)/2

data = readfits(arc,/silent)
nf0 = n_elements(data[*,0])

readcol,mask,silent=1,format='i,i',row,fib
nf1 = n_elements(fib)

index = row - 1

readcol,ptow,silent=1,format='d,d,d,d,d',c0,c1,c2,c3,c4
if (n_elements(c4) ne nf1) then stop

wavearr = dblarr(nf0,nf1)
for jj=0,nf1-1 do begin
    for k=0,nf0-1 do begin
        wavearr[k,jj] = c0[jj]+c1[jj]*k+c2[jj]*k*k+c3[jj]*k*k*k+c4[jj]*k*k*k*k
    endfor
endfor
disp = median(c1)
wi = median(c0)
wave = dblarr(2048)
for jj=0,2047 do wave[jj] = wi + (jj * disp)

nf2 = n_elements(wavelist)
wave0 = wavelist-10
wave2 = wavelist+10

values = fltarr(nf2,nf1,5) ; wavelength,fiber#,(wave,deltawave,line_intensity,FWHM,km/sec)

print,'Fitting gaussians to '+arc
if (splot eq 'on') then window,0,retain=2
for jj=0,nf1-1 do begin ;a loop through fibers
    onew = wavearr[*,jj]
    if (total(onew) eq 0) then goto, jump1
    fiber = total(data[*,index[jj]-step:index[jj]+step],2)
    for kk=0,nf2-1 do begin ;a loop through each wavelength found by PeakFinder
        i = where(onew gt wave0[kk] and onew lt wave2[kk],count)
        if (count gt 9) then begin ;avoids lines at the edge of the chip
            pd = fiber[i]
            pw = onew[i]
            temp = gaussfit(pw,pd,A,nterms=6)
            if (splot eq 'on' and jj eq 39) then begin
                plot,pw,pd,thick=2,title='Fiber: '+strn(jj+1)+' Line: '+strn(kk+1)
                oplot,pw,temp,linestyle=3,thick=2
                pause
            endif
            if (splot eq 'on' and jj eq 209) then begin
                plot,pw,pd,thick=2,title='Fiber: '+strn(jj+1)+' Line: '+strn(kk+1)
                oplot,pw,temp,linestyle=3,thick=2
                pause
            endif
;            if (splot eq 'on' and jj gt 0) then begin
;                plot,pw,pd,thick=2,title='Fiber: '+strn(jj+1)+' Line: '+strn(kk+1)
;                oplot,pw,temp,linestyle=3,thick=2
;                wait,0.1
;            endif
            fwhmA = 2*sqrt(2*alog(2))*A[2]
            values[kk,jj,0] = A[1]
            values[kk,jj,1] = A[1]-wavelist[kk]
;            values[kk,jj,2] = weights[kk]
            values[kk,jj,3] = fwhmA
            values[kk,jj,4] = 299792.458*(A[2]*disp/A[1])
        endif
    endfor
    jump1: ;skips over dead fibers
endfor

i = where(values[0,*,0] eq 0)
if (i[0] ne -1) then values[*,i,*] = -1.0

out = values

return,out

END
