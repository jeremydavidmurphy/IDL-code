pro fwhm

readcol, 'd1.dat', f, w, t, t, t, fwhm

sel1=where(abs(w-3983.9388) le 1.0)
sel3=where(abs(w-4678.1474) le 1.0)
sel4=where(abs(w-4799.9080) le 1.0)
sel5=where(abs(w-4916.0661) le 1.0)
sel6=where(abs(w-5085.8173) le 1.0)
sel7=where(abs(w-5154.6616) le 1.0)
sel8=where(abs(w-5769.6000) le 1.0)

lambda=[median(w[sel1]),median(w[sel3]),median(w[sel4]),median(w[sel5]),median(w[sel6]),median(w[sel7]),median(w[sel8])]
mfwhm=[median(fwhm[sel1]),median(fwhm[sel3]),median(fwhm[sel4]),median(fwhm[sel5]),median(fwhm[sel6]),median(fwhm[sel7]),median(fwhm[sel8])]

window, 0, retain=2, xsize=500, ysize=500

entry_device = !d.name
set_plot , 'ps'
device, /encapsul, filename='virusp_fwhm.ps', xsize=20, ysize=20, /color


plotsym, 0, 2.0, /fill
plot, lambda, mfwhm, psym=8, yrange=[3.5, 7.0], xtitle='Wavelength (Angstroms)',$
  ytitle='FWHM (Angstroms)', charsize=1.5

fit=linfit(lambda, mfwhm)
oplot, findgen(8000), fit[0]+fit[1]*findgen(8000), linestyle=1
oplot, w[sel1], fwhm[sel1], psym=3
oplot, w[sel3], fwhm[sel3], psym=3
oplot, w[sel4], fwhm[sel4], psym=3
oplot, w[sel5], fwhm[sel5], psym=3
oplot, w[sel6], fwhm[sel6], psym=3
oplot, w[sel7], fwhm[sel7], psym=3
oplot, w[sel8], fwhm[sel8], psym=3

device, /close_file
set_plot, entry_device


stop
end

