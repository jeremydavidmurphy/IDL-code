;This routine is used to sum up all the fibers (or a region of fibers)
;for a given frame. The output is a 1-D spectra and the intention is
;to use these frames to determine the best sky-subtraction for a given frame.

;A section of the data can be used, rather than the entire frame. This
;is set by the F1 and FN parameters. THESE ARE NOT THE ACTUAL FIBER
;NUMBERS BUT RATHER THE ROWS FROM THE BOTTOM OF THE DATA FRAME (ref: ds9 window).
PRO sumframe, ForF, wi, wd, ptow, OUTNAME=outname, FRAME=file, F1=f1, FN=fn

n0 = n_elements(file)
if (n0 eq 0) then begin
    file = ''
    print,'Enter the name of the file (w/o the fits):'
    read,file
endif
file = file+'.fits'
n0 = n_elements(outname)
if (n0 eq 0) then begin
    outname = ''
    print,'Enter the output name (w/o the fits):'
    read,outname
endif
outname = outname+'.fits'

data = readfits(file)
wgts = dblarr(n_elements(data[0,*]))
for j=0,n_elements(data[0,*])-1 do begin
    wgts[j] = median(data[*,j])
endfor
m = max(wgts)
wgts = wgts/m

adata = realign(data,ptow,wi,wd,ForF)

n1 = n_elements(f1)
nn = n_elements(fn)
if(n1 ne 0) or (nn ne 0) then begin
    if (n1 eq 0) and (nn ne 0) then adata = adata[*,0:fn-1]
    if (n1 ne 0) and (nn eq 0) then adata = adata[*,f1-1:*]
    if (n1 ne 0) and (nn ne 0) then adata = adata[*,f1-1:fn]
endif

outdata = biweight(adata,NORM=wgts)

writefits,outname,outdata

stop
end
    
