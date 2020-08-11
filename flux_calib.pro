PRO flux_calib, datafile

airmass = 1.1
wi = 3585.0
wdisp = 1.125

spec = readfits(datafile,h)
n0 = n_elements(spec)

window,2,retain=2

readcol,'/scratch/grad78/murphy/research/m87/work/m87f/fluxfactor_-8.dat',f='x,f,f,f,x,x',w,c,x

wave = fltarr(n0)

for j=0,n0-1 do wave[j] = wi + (j * wdisp)

corr = interpol(c,w,wave)
xcorr = interpol(x,w,wave)

spec1 = spec * corr

trans = fltarr(n0)

for j=0,n0-1 do trans[j] = 10^(0.4 * xcorr[j] * airmass)

spec2 = spec1 * trans

plot,wave,spec,xrange=[3600,5500],/xstyle,/ynozero,$
  title='Uncorrected Spectra'
pause

plot,wave,spec1,xrange=[3600,5500],/xstyle,/ynozero,$
  title='Spectra corrected for VIRUS-P'
pause

plot,wave,spec2,xrange=[3600,5500],/xstyle,/ynozero,$
  title='Spectra corrected for atmospheric transmission'
pause


stop
END
