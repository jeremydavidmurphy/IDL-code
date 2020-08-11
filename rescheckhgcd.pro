; This version is specific for the calcRES.pro routine. It does the
; fitting of the gaussians to each emission line in each fiber. The
; plotting turned on then generates the visual inspection of the
; quality of fits.

; The primary difference between this and rescheck2F is that in
; rescheck2F, the routine runs a peak-finding routine to locate the
; lines. As the arcs are well known, this routine doesn't do the
; search, but just looks for lines within +/- 10 A of the lines in the
; list 'wavelist'

FUNCTION rescheckhgcd, arc,ptow,mask,nwave,Splot=splot

if (nwave eq 10) then $
wavelist = [4046.5469,4077.8403,4358.3262,4678.1558,4799.9038,$
            4916.0962,5085.8110,5460.7397,5769.5972,5790.6348]
if (nwave eq 9) then $
wavelist = [4046.5469,4077.8403,4358.3262,4678.1558,4799.9038,$
            4916.0962,5085.8110,5460.7397,5769.5972];this was used for feb08, where this line runs off the chip.
if (nwave eq 8) then $
wavelist = [4046.5469,4077.8403,4358.3262,4678.1558,4799.9038,$
            5085.8110,5460.7397,5769.5972];oct07
if (n_elements(Splot) eq 0) then splot = 'off'

nprofile = 5 ;the number of rows in a fiber profile
step = (nprofile-1)/2

data = readfits(arc,/silent)
nf0 = n_elements(data[*,0])
nf3 = n_elements(data[0,*])

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
wave = dblarr(nf3)
for jj=0,nf3-1 do wave[jj] = wi + (jj * disp)

nf2 = n_elements(wavelist)
wave0 = wavelist-10
wave2 = wavelist+10

values = fltarr(nf2,nf1,5) ; wavelength,fiber#,(wave,deltawave,line_intensity,FWHM,km/sec)

print,'Fitting gaussians to '+arc+'...'
pcntr = 9
if (splot eq 'on') then begin
    window,0,retain=2
    device,decomposed=0
    loadct,0,/silent
endif
for jj=0,nf1-1 do begin ;a loop through fibers
    onew = wavearr[*,jj]
    if (total(onew) eq 0) then begin
        print,'Jumping dead fiber #'+strn(jj+1)+'...'
        goto, jump1             ;a step over dead fibers
    endif
    fiber = total(data[*,index[jj]-step:index[jj]+step],2)
    for kk=0,nf2-1 do begin ;a loop through each wavelength found in the list
        i = where(onew gt wave0[kk] and onew lt wave2[kk],count)
        if (count gt 8) then begin ;avoids lines at the edge of the chip
            pd = fiber[i]
            pw = onew[i]
            temp = gaussfit(pw,pd,A,nterms=3)
;            if (splot eq 'on' and jj eq 200) then begin
;                plot,pw,pd,thick=2,title='Fiber: '+strn(jj+1)+' Line: '+strn(kk+1)
;                oplot,pw,temp,linestyle=3,thick=2
;            endif
;            if (splot eq 'on' and jj eq 245) then begin
;                plot,pw,pd,thick=2,title='Fiber: '+strn(jj+1)+' Line: '+strn(kk+1)
;                oplot,pw,temp,linestyle=3,thick=2
;                pause
;            endif
            if (splot eq 'on' and jj eq pcntr) then begin
                plot,pw,pd,thick=2,title='Fiber: '+strn(jj+1)+' Line: '+strn(wavelist[kk])
                oplot,pw,temp,linestyle=3,thick=2
                wait,0.1
                if (kk eq nf2-1) then pcntr = pcntr + 10
            endif
            fwhmA = 2*sqrt(2*alog(2))*A[2]
            if fwhma lt 2 then stop
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
