pro collapse_sky, skylist

;This collapses the name_pefy.fits files into 1-D spectra.
;It accepts a list of frames, then renames the collapsed data as
;name_skyC.fits

; MODIFIED (Aug 3rd, 2010): Added in the wavelength shift (via
; realignF and inclusion of the ptow.dat file. The SKYLIST is now in
; the form:

; jm####_pefy.fits   ptow_month_night.dat   mask_month_night.dat

;***********************************************************
wi = 3550.0
wdisp = 1.125
factor = 3.0
app = 5
plotting = 'on'
;***********************************************************

readcol,skylist,format='a,a,a',list,ptow,mask

n1 = n_elements(list)

t = readfits(list[0])
n2 = n_elements(t[*,0])
n3 = n_elements(t[0,*])
skyplot.ps
wave = fltarr(n2*factor)
wave = findgen(n2*factor)*(wdisp/factor) + wi

if (plotting eq 'on') then window,2,retain=2

for j=0,n1-1 do begin
    onesky = readfits(list[j],h)
    print,'Working on sky frame: '+list[j]
    wzp = sxpar(h,'WAVEZP',count=cnt)
    if (cnt eq 0) then wzp = 0.0
    helioc = sxpar(h,'VHELIO',count=cnt)
    if (cnt eq 0) then helioc = 0.0
    asky = realignf(onesky,ptow[j],mask[j],wi,wdisp,app,factor=factor,wzp=wzp,helioc=helioc)
    msky = median(asky,dim=2)
    if (plotting eq 'on') then begin
        loadct,0,/silent
        plot,wave,msky,/ynozero,title='Collapsed Sky Spectra: '+list[j],$
          yrange=[0,median(msky)+15],/ys
    endif
    if (j eq 0) then begin
        w1out = wdisp/factor
        headout = h
        sxaddpar,headout,'CRVAL1',wi,' Initial wavelength'
        sxaddpar,headout,'CDELT1',w1out,' Wavelength dispersion term'
        sxaddpar,headout,'COMMENT','  This spectra has been interpolated to a common wavelength solution'
    endif
    o = strsplit(list[j],'_',/extract)
    outname = o[0]+'_skyC.fits'
    out = msky
    writefits,outname,out,headout
endfor

wait,3
wdelete

stop
end
