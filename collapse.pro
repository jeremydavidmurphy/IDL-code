; The code collapses VIRUS-P data. It allows for a linear
; interpolation via the realigncF.pro fuction. The user will be
; prompted for this step.

pro collapse, dataname, maskname, type, PTOW=ptowname, $
              ONE=one, WZP=zero, WI=wi, DISP=disp

; this function collapses VIRUS-P data. It requires the mask file. The
; data files need to have an 'e' in the suffix (meaning they've been extracted).

; TYPE: This is a string of either "total", "median", or "mean" and
; dictates how the fiber is collapsed.

; PTOW: The name of the ptow file. If this is included then the output
; file will be aligned to a linear wavelength solution.

; ONE: If the keyword ONE is used (you can put any string in here),
; the file returned is a combined set of a range of fibers. You will
; be prompted for which fibers later in the routine. It is output as
; "oneC.fits".

; WZP: If the "wzp" callword is set to 'zero' then the header will be
; searched for the WAVEZP parameter. This will zero off the 5577
; oxygen sky line.

; WI: The initial wavelength you want to linearize your data to. You
; will be prompted for this, if not included and you are aligning.

; DISP: The dispersion of your data. You will be prompted for this one
; as well.

;*************************************************************
aperture = 5.0 ;the size of the extraction aperture
step = (aperture-1)/2.0
;*************************************************************

data = readfits(dataname,header,/silent)

readcol,maskname,silent=1,format='i,x',index
index = index-1

n0 = n_elements(index)
n1 = n_elements(data[*,0])

dataout = fltarr(n1,n0)

badf = where(data eq -666,count)
if (count gt 0) then begin
    data[where(data eq -666)] = 0.0
    print,'You have -666 flags in your data!'
    print,'These will be set to zero...'
endif

;The collapse is made...
for j=0,n0-1 do begin ;a loop through each fiber
    row = data[*,index[j]-step:index[j]+step]
    if (type eq 'total') then dataout[*,j] = total(row,2)
    if (type eq 'median') then dataout[*,j] = median(row,dimension=2,/even)
    if (type eq 'mean') then begin
        for k=0,n1-1 do dataout[k,j] = mean(row[k,*])
    endif
endfor

if (n_elements(ptowname) eq 0) then begin ;the file is just written out and you are done
    temp = strsplit(dataname,'.',/extract)
    temp = temp[0]
    nameout = temp+'C.fits'
    print,'The file has been written out as '+nameout
    headout = header
    sxaddpar,headout,'COMMENT','  This spectra was collapsed using collapse.pro'
    writefits,nameout,dataout,headout
    goto,jumpend
endif else begin ;the file is realigned
    
    if (n_elements(wi) eq 0) then begin
        wi = 1.0
        print,'Enter an initial wavelength:'
        read,wi
    endif
    if (n_elements(disp) eq 0) then begin
        wdisp = 1.0
        print,'Enter a dispersion:'
        read,wdisp
    endif
    if (n_elements(wzp) ne 0) then begin
        wzp = sxpar(header,'WAVEZP',count=c)
        if (c ne 1) then wzp = 0.0
        wzp = double(wzp)
    endif
    aligned = realigncf(dataout,ptowname,wi,wdisp,WZP=wzp)
    temp = strsplit(dataname,'.',/extract)
    temp = temp[0]
    nameout = temp+'AC.fits'
    print,'The file has been written out as '+nameout
    headout = header
    swi = string(wi)
    swdisp = string(wdisp)
    sxaddpar,headout,'COMMENT','  This spectra was collapsed using collapse.pro'
    sxaddpar,headout,'COMMENT','  This spectra has been interpolated to a common wavelength solution'
    sxaddpar,headout,'CRVAL1',wi,' Initial wavelength',after='IMAGETYP'
    sxaddpar,headout,'CDELT1',wdisp,' Wavelength dispersion term',after='CRVAL1'
    writefits,nameout,aligned,headout
endelse
    
if (n_elements(one) ne 0) then begin
    fiber1 = 1
    fiber2 = 2
    print,'Enter a lower fiber number:'
    read,fiber1
    fiber1 = fiber1 - 1
    print,'Enter an upper fiber number:'
    read,fiber2
    fiber2 = fiber2 - 1
    dif = fiber2-fiber1+1
    fibers = indgen(dif)+fiber1
    print,'The fibers going into the combination are:'
    print,fibers+1
    piece = aligned[*,fibers]
    if (type eq 'total') then pieceout = total(piece,2)
    if (type eq 'median') then pieceout = median(piece,dimension=2)
    if (type eq 'mean') then begin
        pieceout = fltarr(n1)
        for k=0,n1-1 do pieceout[k] = mean(piece[k,*])
    endif
    headout = header
    sxaddpar,headout,'CRVAL1',wi,' Initial wavelength'
    sxaddpar,headout,'CDELT1',wdispF,' Wavelength dispersion term'
    sxaddpar,headout,'COMMENT','  This spectra was collapsed using collapse.pro'
    sxaddpar,headout,'COMMENT','  This spectra has been interpolated to a common wavelength solution'
    sxaddpar,headout,'COMMENT','  This spectra is a combination of a range of fibers'
    writefits,'oneC.fits',pieceout,headout
endif

jumpend:

stop
END
