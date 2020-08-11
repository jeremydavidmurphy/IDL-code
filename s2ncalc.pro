; This routine is used to calculate the S/N of your spectra based on
; the output of pipe2.pro. This value is required in the monte carlo
; simulations used to estimate the error of your losvd's. It generates
; the run file for HOKU.

pro s2ncalc, wi, wd, RUN=runfile

; WI and WD are the initial wavelength and wavelength dispersion of
; your data. This should be the same value that you put into PIPE2.PRO.

; This code requires a few files to exist in the calling directory
; s2n.list: This is a list of all the s2n.fits files created by
; PIPE2.PRO
; region.list: This is a list of the wavelength (matching that of all
; the run files) for each bin. Example below:
; GB   4195  4585
; HB   4455  4945
; MGw  5000  5410
; MGwo 5000  5410
; FE   5280  5790

; The 'run3' file (or the general run file) is just the run file you
; used to generate the pfitlov output (i.e. rfitlov).

scale = 0.3 ; a factor that gets multiplied by your s2n calc, as it seems too noisy...

readcol,silent=1,'s2n.list',format='a',fitsfiles ; the list of name_s2n.fits files
readcol,silent=1,'region.list',format='a,i,i',regions,wr1,wr2

if (n_elements(runfile) eq 0) then begin
    runfile = ''
    print,'Enter the name of your run file...'
    read,runfile
endif
readcol,silent=1,runfile,format='x,a,i,i,f,i,i,a',b1,w1,w2,z,s,sm,reg

test = readfits(fitsfiles[0])
n0 = n_elements(fitsfiles)
n1 = n_elements(test)

wavedata = fltarr(n1)
for j=0,n1-1 do wavedata[j] = wi + (wd *j)

bins = strarr(n1)

n2 = n_elements(b1)
r2t = strarr(n2)
rms = fltarr(n2)

for j=0,n0-1 do begin
    file = fitsfiles[j]
    bin = strsplit(file,'_',/extract)
    bin = bin[0]
    s2n = readfits(file,/silent)
    ibin = where(b1 eq bin)
    if (ibin[0] ne -1) then begin
        for k=0,n_elements(ibin)-1 do begin
            wave1 = w1[ibin[k]]
            wave2 = w2[ibin[k]]
            iw1 = where(wavedata le wave1)
            iw2 = where(wavedata gt wave2)
            iw1 = iw1[n_elements(iw1)-1]
            iw2 = iw2[0]
            piece = s2n[iw1:iw2]
            rms[ibin[k]] = float(1/mean(piece)) * scale
;**************************************************
; a total fudge to reign in the low and high values...

 ;           if (rms[ibin[k]] lt 0.09) then rms[ibin[k]] = 0.5 * rms[ibin[k]] $
 ;           else rms[ibin[k]] = 0.06
;            if (rms[ibin[k]] lt 0.007) then rms[ibin[k]] = rms[ibin[k]] * 2.0
;            if (rms[ibin[k]] lt 0.007) then rms[ibin[k]] = 0.007
            if (rms[ibin[k]] gt 0.09) then rms[ibin[k]] = 0.07
            if (rms[ibin[k]] lt 0.00) then rms[ibin[k]] = 0.08
            if (rms[ibin[k]] lt 0.004) then rms[ibin[k]] = 0.004
;**************************************************
        endfor
    endif
endfor

free_lun,5
openw,5,'run.out'
for j=0,n2-1 do begin
    print,b1[j],w1[j],w2[j],z[j],s[j],sm[j],reg[j],rms[j]
    printf,5,'rmcfit  ',b1[j],w1[j],w2[j],z[j],s[j],sm[j],reg[j],rms[j],$
  format='(a8,a7,2x,i4,2x,i4,2x,f8.6,2x,i3,1x,i2,1x,a4,1x,f9.6)'
endfor
free_lun,5

stop
end
