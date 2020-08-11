; This code reads in an arc,ptow and mask file, and returns the IR for
; that night, based on Gaussian fits to the arc lines.

; COMPILE: Make_Set
; COMPILE: PeakFinder

FUNCTION rescheck2F, arc,ptow,mask,num,Splot=splot

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
;and the center of the line, as measu

if (n_elements(Splot) eq 0) then splot = 'off' else splot='on'

nprofile = 5 ;the number of rows in a fiber profile
step = (nprofile-1)/2

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

wavexplore = fltarr(num,2,nf1)

for jj=0,nf1-1 do begin
    one = total(data[*,index[jj]-step:index[jj]+step],2)
    find = PeakFinder(one,wave[*,jj])
    done = 'go'
    cutoff = 0.00001
    repeat begin
        ilim = where(find[4,*] gt cutoff) ;indexes those values above the cutoff
        if (n_elements(ilim) gt num) then cutoff = cutoff + 0.00001 else $
          done = 'stop'
    endrep until (done eq 'stop')
    wavexplore[*,0,jj] = find[1,ilim]
    wavexplore[*,1,jj] = find[4,ilim]
endfor

strength = median(wavexplore[*,1,*],dim=3)
wavelist = median(wavexplore[*,0,*],dim=3)
nf2 = n_elements(wavelist)
wave0 = wavelist-12
wave2 = wavelist+12

collapsed = fltarr(nf0,nf1)
values = fltarr(nf2,nf1,5) ; wavelength,fiber#,(wave,deltawave,line_intensity,FWHM,km/sec)

print,'Fitting gaussians to '+arc
if (splot eq 'on') then window,2,retain=2
for jj=0,nf1-1 do begin ;a loop through fibers
    onew = wave[*,jj]
    if (total(onew) eq 0) then goto, jump1
    fiber = total(data[*,index[jj]-step:index[jj]+step],2)
    for kk=0,nf2-1 do begin ;a loop through each wavelength in the arcs
        i = where(onew gt wave0[kk] and onew lt wave2[kk],count)
        if (count gt 9) then begin ;avoids lines at the edge of the chip
            pd = fiber[i]
            pw = onew[i]
            temp = gaussfit(pw,pd,A,nterms=4)
            if (splot eq 'on' and jj eq 0) then begin
                plot,pw,pd,thick=2,title=strn(kk+1)
                oplot,pw,temp,linestyle=3,thick=2
                pause
            endif
            if (splot eq 'on' and jj gt 0) then begin
                plot,pw,pd,thick=2,title='Fiber: '+strn(jj+1)+'Line: '+strn(kk+1)
                oplot,pw,temp,linestyle=3,thick=2
                wait,0.1
            endif
            fwhmA = 2*sqrt(2*alog(2))*A[2]
            values[kk,jj,0] = A[1]
            values[kk,jj,1] = A[1]-wavelist[kk]
            values[kk,jj,2] = strength[kk]
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
