
pro collapse_all, framelist

; compile realignf.pro

;This collapses an entire VIRUS-P frame into 1-D spectra.
;It accepts a list of frames, then renames the collapsed data as
;nameC.fits. 

; THIS CODE DOES NOTHING FANCY. IT COLLAPSES EACH FIBER BY TAKING THE
; TOTAL OVER THE 5 PIXEL APERTURE. IT THEN TAKES THE MEDIAN OR MEAN
; OVER ALL THE FIBERS, DEPENDING ON THE KEYWORD YOU SET BELOW.

; FRAMELIST: A list of 3 columns that gives the fits file name, the
; ptow and mask file for that frame.

; MODIFIED (Aug 3rd, 2010): Added in the wavelength shift (via
; realignf and inclusion of the ptow.dat file.

;***********************************************************
wi = 3550.0
wdisp = 1.125
app = 5
plotting = 'on'
morm = 'median' ;change to 'mean' if you want to take the average over all fibers. (this has the drawback of getting hit with cosmic rays and not recovering very well...)
;***********************************************************

readcol,framelist,format='a,a,a',frames,ptow,mask

n1 = n_elements(frames)

t = readfits(frames[0])
n2 = n_elements(t[*,0])
n3 = n_elements(t[0,*])

wave = findgen(n2)*wdisp + wi

if (plotting eq 'on') then window,2,retain=2

openw,5,'collapsed_values.txt'

for j=0,n1-1 do begin
    frame = readfits(frames[j],h)
    print,'Working on frame: '+frames[j]
    cntr = 0
    out = 'lost'
    repeat begin
        hh = strsplit(h[cntr],' ',/extract)
        if (hh[0] eq 'WAVEZP') then out = 'found' else cntr = cntr + 1
    endrep until (out eq 'found')
    wzp = float(hh[2])
    aframe = realignf(frame,ptow[j],mask[j],wi,wdisp,app,factor=1,wzp=wzp,helioc=0.0)

    ; now the frame is collapsed to single spectra for each fiber
    readcol,mask[j],format='i,x',y,silent=1
    i = y - 1.0
    nm = n_elements(i)
    col = fltarr(n2,nm)
    for k=0,nm-1 do begin
        for kk=0,n2-1 do col[kk,k] = total(aframe[kk,i[k]-2:i[k]+2])
    endfor

    if (morm eq 'median') then mframe = median(col,dim=2,/even)
    if (morm eq 'mean') then begin
        mframe = fltarr(n2)
        for k=0,n2-1 do mframe[k] = mean(col[k,*])
    endif
    if (plotting eq 'on') then begin
        yup = median(mframe[0:1500])*2.0
        ydn = 0.0
        plot,wave,mframe,/ynozero,title='Collapsed Spectra: '+frames[j],$
          yrange=[ydn,yup]
    endif
    o = strsplit(frames[j],'.',/extract)
    outname = o[0]+'C.fits'
    out = [[wave],[mframe]]
    writefits,outname,out
    o1 = string(total(mframe))
    o2 = string(median(mframe,/even))
    o3 = string(median(mframe[1000:1400],/even))
    o4 = string(mean(mframe))
    printf,5,[o[0],o1,o2,o3,o4]
endfor
free_lun,5

stop
end
