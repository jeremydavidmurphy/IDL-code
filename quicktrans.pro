PRO quicktrans

; This code is for a quick reduction of transmission data for VIRUS
; optical fibers. It works from 5 lists. They can be of any length
; (the length is the number of wavelengths) but must have the
; following names:
; bllist: the list of baseline data frames
; bldlist: the list of baseline dark frames
; flist: the list of fiber data frames
; fdlist: the list of fiber dark frames
; wlist: the wavelengths

readcol,'bllist',silent=1,format='a',blfiles
readcol,'bldlist',silent=1,format='a',bldfiles
readcol,'flist',silent=1,format='a',ffiles
readcol,'fdlist',silent=1,format='a',fdfiles
readcol,'wlist',silent=1,format='i',wave

n0 = n_elements(wave)
trans = dblarr(n0)

for j=0,n0-1 do begin
    print,'Working on wavelength '+strn(wave[j])+'...'
    bl = readfits(blfiles[j],/silent)
    bld = readfits(bldfiles[j],/silent)
    fib = readfits(ffiles[j],/silent)
    fibd = readfits(fdfiles[j],/silent)
    bls = bl - bld
    fibs = fib - fibd
    trans[j] = total(fibs) / total(bls)
endfor

window,0,retain=2
plot,wave,trans,psym=-1,title='Quick Transmission Plot',charsize=1.2,$
  symsize=1.5,xtitle='Wavelength (A)',ytitle='Transmission'

stop
END
