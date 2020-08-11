; This code reads in an arc,ptow and mask file, and returns the IR for
; that night, based on Gaussian fits to the arc lines.

; The difference between 3F and 2F is that in this code the peak
; finding is found is a more robust way (aligning, then collapsing the
; arc frame).

; COMPILE: Make_Set
; COMPILE: PeakFinder
; COMPILE: realignF

FUNCTION rescheck3F, arc,ptow,mask,num,Splot=splot

; Modified on Sept 2010: The difference between this version and
; rescheckF is that now the code PeakFinder is used to locate all the
; available lines. It does this by searching fiber #99 for lines and
; then using this information (i.e. the location of the peaks for this
; fiber) for all the other fibers. The original version used an
; internal list of lines.

;NUM: The number of lines to find and fit with arclines.

;NOTE: The output is a 3D array of the following form: ARRAY[X,Y,Z]
;ARRAY[X]: wavelength. This will be of the same dimension as NUM
;ARRAY[Y]: fiber number
;ARRAY[Z]: 4 values: wavelength, delta_wave, FWHM, km/sec
;Delta_wave is the difference between the center of the Gaussian fit
;and the center of the line.

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

; This first bit locates the peaks...

aligned = realignF(data,ptow,mask,wi,disp,nprofile)
alignedC = dblarr(n_elements(aligned[*,0]),nf1)

for jj=0,nf1-1 do begin
    alignedC[*,jj] = total(aligned[*,index[jj]-step:index[jj]+step],2)
endfor

ONE = median(alignedC,dim=2,/even)

find = PeakFinder(ONE,wave)
done = 'go'
cutoff = 0.00001
repeat begin
    ilim = where(find[4,*] gt cutoff) ;indexes those values above the cutoff
    if (n_elements(ilim) gt num) then cutoff = cutoff + 0.00001 else $
      done = 'stop'
endrep until (done eq 'stop')

wavelist = find[1,ilim]
weights = find[4,ilim]
 
nf2 = n_elements(wavelist)
wave0 = wavelist-10
wave2 = wavelist+10

values = fltarr(nf2,nf1,5) ; wavelength,fiber#,(wave,deltawave,line_intensity,FWHM,km/sec)

print,'Fitting gaussians to '+arc
if (splot eq 'on') then window,2,retain=2
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
            values[kk,jj,2] = weights[kk]
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
