;This routine is used to sum up all the fibers (or a region of fibers)
;for a given frame. The output is a 1-D spectra and the intention is
;to use these frames to determine the best sky-subtraction for a given frame.

;A section of the data can be used, rather than the entire frame. This
;is set by the F1 and FN parameters. THESE ARE NOT THE ACTUAL FIBER
;NUMBERS BUT RATHER THE ROWS FROM THE BOTTOM OF THE DATA FRAME (ref: ds9 window).

;PRO sumframe, ForF, wi, wd, ptow, OUTNAME=outname, FRAME=file, F1=f1, FN=fn

;This has been turned into a list driven loop so that many files can
;be done at once. The old format is above. The list has 8 columns.
;EXAMPLE:
;FILE                  OUTNAME          WI      WD    PTOW         ForF     F1     FN
;jm0077aa_pefsmc.fits  jm0077aaTT.fits  3530.0  1.12 ptow_jan08_n1.dat  force  201  247

;--------------------------------------------------------------
PRO sumframe_v2, list
;--------------------------------------------------------------

;When run with the biweight over all fibers, this takes ~10s per frame.

ans=''
print,'Plot the output? ("y" or "n")'
read,ans
if(ans eq 'n') then begin
    print,''
    print,'Plotting is OFF!'
    print,''
    pause
endif

readcol,list,format='a,a,f,f,a,a,i,i',file,outname,wi,wd,ptow,forf,f1,fn
n0 = n_elements(file)
test = readfits(file[0])
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*])

if (ans eq 'y') then begin
    window,1,retain=2,xsize=900
    loadct,0
endif

for k=0,n0-1 do begin
    data = readfits(file[k])
    if(n_elements(data) eq 1) then begin
        print,'File '+file[k]+' was not found!'
        print,'Skipping to the next file in the list.'
        goto,jumpend
    endif
    
    adata = realign(data,ptow[k],wi[k],wd[k],ForF[k])

    ff1 = f1[k]
    ffn = fn[k]
    n1 = n_elements(ff1)
    nn = n_elements(ffn)
    if(n1 ne 0) or (nn ne 0) then begin
        if (n1 eq 0) and (nn ne 0) then adata = adata[*,0:ffn-1]
        if (n1 ne 0) and (nn eq 0) then adata = adata[*,ff1-1:*]
        if (n1 ne 0) and (nn ne 0) then adata = adata[*,ff1-1:ffn-1]
    endif

    n3 = n_elements(adata[0,*])
    wgts = dblarr(n3)
    for j=0,n3-1 do begin
        wgts[j] = median(adata[*,j])
    endfor
    m = max(wgts)
    wgts = wgts/m

    outdata = biweight(adata,NORM=wgts)
;    for j=0,n3-1 do adata[*,j] = adata[*,j]*wgts[j]
;    outdata = total(adata,2)

    writefits,outname[k],outdata
    yup = max(outdata[20:2020])
    if(ans eq 'y') then begin
        plot,outdata,title='THE OUTPUT FOR '+file[k],charsize=1.5,$
          yrange=[0,yup],xrange=[0,2048],xstyle=1
    endif

    jumpend:

endfor

end
    
