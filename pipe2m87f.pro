; The second version of pipe2 has a few modifications. It allows for
; either the biweight, or a weighted average combination with sigma
; clipping. It also allows for the instrumental response and/or
; transmission loss due to airmass effects to be backed out of the
; data. Finally, it allows the wavelength solution to be
; super-sampled. This mitigates the smoothing effects that occur with
; linear interpolation.

; This is the second half of the pipeline and assumes you've run
; PIPE1.PRO in the same directory you're running PIPE2.PRO
; from. This is necessary as PIPE1.PRO generates the name.R files
; required by PIPE2.PRO. This procedure works from lists with the
; ultimate goal of returning the final, aligned, and combined file for
; a given pointing. Pipe2 is a reworking of alignall.pro, realign.pro
; and the combiwgt.pro procedures. 

; pipe2 returns air wavelengths. to switch to vacuum, set the VAC
; keyword equal to 'yes' if you want the wavelength solution to be
; converted to vacuum wavelengths.
;*****************************************************************
pro pipe2m87f, LIST=listoflists, BIN=yourbinlist,$
              COMB=COMB, IR=IR, TRANS=TRANS,$
              SKIPINTERPOL=skipintpol, VAC=vac
;*****************************************************************

; THE ONLY DIFFERENCE HERE IS THAT THERE ARE SO MANY FILES I CAN'T
; ALLOCATE ENOUGH MEMORY. SO, THE FILES WILL BE COLLAPSED INTO SINGLE FI
;COMPILE :realignf,biweightf,meanwclipf

; EXAMPLE OF CALLING SEQUENCE:
; pipe2, list='N4889data.list', COMB='biwt', IR='yes', TRANS='no'
; In general, if you've named your list of bins, 'bin.list' and your
; data list, 'data.list', this thing will run with just the call to pipe2_v2.

; MODIFIED ON JULY 17, 2011: The code was modified to handle very
; large arrays (e.g. M87f).

; set VAC equal to anything and the wavelength solution will be
; converted to vacuum wavelengths. This is done within the routine
; realignF, so see that piece of code for details.

;***************************************************************************
wi = 3550.0 ;(M87:3545, M87f:3580, N3842:3575, N4472:3550, N2832:3540) ;The initial wavelength and dispersion values for the linear interpolation.
wdisp = 1.125 ;(M87:1.125) NOTE: The 'factor' term makes this a 1.125/3.0 = 0.375 dispersion term
factor = 3.0 ;this allows for super-sampling of the wavelength solution by changing the 'wdisp' value to wdisp/factor
aperture = 5 ;The size of the fiber profile, in pixels.
lastplot = 'on' ;This just plots the final spectra for each bin once complete.
filesuffix = 'pefsm' ;This dictates whether you're combining the pefs or pefsm files
pointarr = ['a1','a2','a3','a','b','c','d','e','f','g','h','d1','d2','d3'] ;all the available pointings
focalreducer = 'new'
;focalreducer = 'old'
;IFU = 'vp1'
IFU = 'vp2'
COMP = 'grad78' ;change to grad78 when running on your grad box (dictates where to find the flux correction files)
tsci = 1800.0 ;***
tsky = 300.0 ;***
tfibers = 229.0 ;***
;***************************************************************************

;The listoflists is just that. A list of different data lists. The
;only purpose of this is to allow different mixes of sky-subtraction
;to be run sequentially and without supervision.
;Each of the datalists must have the form below. The listoflists will
;look something like
;  M87cc.list
;  M87bc.list
;  etc.

;where each line is the name of a text file. Each of these text files
;must have the structure shown in the section below

;If the LIST keyword call isn't used, the routine searches for a
;list called "data.list" 

;The datalist:
;THE FORMAT FOR THE DATALIST (i.e. M87cc.list)
;The list is composed of 4 columns. The first is the data
;file, sans any pefsm.fits. The second indicates the pointing,
;corresponding to the binlist letters. The third column is the
;corresponding ptow.dat AND mask.dat files. This will be the name you
;gave VACCINE which gets used in the naming convention for the data
;for that run. PIPE2.PRO will append the necessary 'mask' and 'ptow'
;as needed. The 4th column are the "name.R" files that are output from 
;PIPE1.PRO. 

;This is a list of the form:
;  jm0123cb_   a   _jan08_n2.dat    jm0123.R   jm0122_   jm0124_
;  jm0125cd_   a   _jan08_n2.dat    jm0125.R   jm0124_   jm0126_
;  etc.

; (NOTE: The lists have been modified to include the two sky nod
; frames.)

;This list must include ALL THE DATA for a given galaxy, INCLUDING ALL THE
;POINTINGS, that you want to go into the final data. This is 
;necessary due to the overlapping pointings.

;YOURBINLIST:
;This is just a list of the bin lists. So, the bin.list should look
;like,
;  bin1601
;  bin1602
;  etc.
;With bin1601 being a text file of the form
;  203_a
;  218_a
;  1_c
;  2_c
;  etc.
; (this routine can handle any number of -1 values in the list, which
; is just an artifact of a step in pipe1)

;IF THE BIN KEYWORD IS NOT USED THE CODE WILL LOOK FOR A FILE NAMED "BIN.LIST"

; COMB: This switch was added on 12/24/2010 to allow the combination
; to be a weighted average. This uses the meanWclipF routine. Accepts 'biwt'
; or 'mean'. If nothing is set, the code will default to biweight.

; IR: ('yes' or 'no') This switch allows you to remove the instrumental response and
; correct for the transmission (based on the calculated airmass). If
; set to 'yes' the IR will be removed. If set to 'no' CCD values will
; be returned with no correction for the IR. If not used, the user
; will be prompted.

; TRANS: ('yes' or 'no') This switch allows you to correct for loss due to
; airmass. THIS SWITCH WILL ONLY WORK IF YOU HAVE THE IR SWITCH ON AS
; WELL.

; SKIPINTERPOL: ('yes' or 'no') The interpolation step is long. Set this keyword to
; 'yes' and the routine will look in the calling directory for the
; name_pefsmA.fits frames. These get created the first time through
; pipe2. If you are re-running pipe2 for some reason, and these files
; got created, you can skip the interpolation step.

;STEPS OF THE REDUCTION:
; All the data and weight arrrays are read in. The weight arrays
; contain values outside of 0 and 1. These are generally data of very
; poor quality, and so the weights are set equal to zero.

; The data is run through realignF.pro, which makes the interpolation
; of the wavelength. This runs on both the data and weight array.

; All the data and weights for a given bin are sorted out. The
; possible pointings are those in the 'pointarr' array, defined
; above. To add a pointing, simply add it to this array.

; MODIFIED ON FEB 8, 2012: The code now handles the errors from the
; sky nods properly. The sky nod frame errors are read in and applied
; during the S/N calculation. This necessitated a change to the format
; of the lists. EVERY PIECE OF THE CODE CHANGED FOR THIS UPGRADE HAS
; BEEN MARKED WITH A *** AT THE END OF THE LINE.

;***************************************************************************

step = (aperture-1)/2
wdispF = wdisp/factor

ans = ''

if (n_elements(vac) gt 0) then vac = 'yes' else vac = 'no'
if (n_elements(IR) eq 0) then begin
    IR = ''
    print,'Would you like the instrumental response'
    print,'removed from the data?'
    print,'ENTER "yes" or "no"'
    read,IR
    if (IR eq 'yes' and n_elements(TRANS) eq 0) then begin
        TRANS = ''
        print,'Do you want to correct for airmass loss?'
        print,'ENTER "yes" or "no"'
        read,TRANS
    endif
endif

if (IR eq 'yes') and (COMP eq 'grad78') then begin
    if focalreducer eq 'old' and IFU eq 'vp1' then $
      readcol,'/home/grad78/murphy/flux_calib/fluxfactor_VP1_OFR_FI1_VPH2.txt',$
      silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
    if focalreducer eq 'old' and IFU eq 'vp2' then $
      readcol,'/home/grad78/murphy/flux_calib/fluxfactor_VP2_OFR_FI1_VPH2.txt',$
      silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
    if focalreducer eq 'new' and IFU eq 'vp2' then $
      readcol,'/home/grad78/murphy/flux_calib/fluxfactor_VP2_NFR_FI1_VPH2.txt',$
      silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
    fluxcor = fluxcor * 1e-17
    
    if (n_elements(TRANS) eq 0) then begin
        TRANS = ''
        print,'Do you want to correct for airmass loss?'
        print,'ENTER "yes" or "no"'
        read,TRANS
    endif
endif

if (IR eq 'yes') and (COMP eq 'danu') then begin
    if focalreducer eq 'old' and IFU eq 'vp1' then $
      readcol,'/home/danu/murphy/research/flux_calib/fluxfactor_VP1_OFR_FI1_VPH2.txt',$
      silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
    if focalreducer eq 'old' and IFU eq 'vp2' then $
      readcol,'/home/danu/murphy/research/flux_calib/fluxfactor_VP2_OFR_FI1_VPH2.txt',$
      silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
    if focalreducer eq 'new' and IFU eq 'vp2' then $
      readcol,'/home/danu/murphy/research/flux_calib/fluxfactor_VP2_NFR_FI1_VPH2.txt',$
      silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
    fluxcor = fluxcor * 1e-17
    
    if (n_elements(TRANS) eq 0) then begin
        TRANS = ''
        print,'Do you want to correct for airmass loss?'
        print,'ENTER "yes" or "no"'
        read,TRANS
    endif
endif

if (n_elements(COMB) eq 0) then COMB = 'biwt'

if (n_elements(listoflists) eq 0) then listoflists = 'data.list'
readcol,silent=1,listoflists,format='A',allists

if (n_elements(skipintpol) eq 0) then skipintpol = 'no'

;The bins are read in
if (n_elements(yourbinlist) eq 0) then readcol,silent=1,'bin.list',format='a',binlist $
else readcol,silent=1,yourbinlist,format='a',binlist
n0 = n_elements(binlist)

for p = 0, n_elements(allists) - 1 do begin ;a loop through each of the lists.
    onelist = allists[p]
    print,''
    print,'Working on list '+onelist+'...'
    print,''
;    readcol,silent=1,onelist,format='a,a,a,a',files,pointing,dat,dotR ;***
    readcol,silent=1,onelist,format='a,a,a,a,a,a',files,pointing,dat,dotR,skyn1,skyn2 ;***

    n3 = n_elements(files) ;the number of files (i.e. number of 20 minute pointings)

    datafilenames = strarr(n3) ;the name of the aligned output files (used for reading in later)
    noisefilenames = strarr(n3)
    skyfilenames = strarr(2,n3) ;***

    if (skipintpol eq 'no') then begin ;the file is realigned and written out

        for j=0,n3-1 do begin ; a loop through each file (the length of the datalist)

            print,''
            print,'Reading in file '+files[j]
            data = readfits(files[j]+filesuffix+'.fits',/silent,dataheader)
            weight = readfits(files[j]+'pefw.fits',/silent,noiseheader)
            sky1 = readfits(skyn1[j]+'pefw.fits',/silent,sky1header) ;***
            sky2 = readfits(skyn2[j]+'pefw.fits',/silent,sky2header) ;***

            izero = where(weight eq 0.0,ct)
            if (ct ne 0) then weight[izero] = 1.0
            noise2 = 1.0 / weight ;THIS IS THE SQUARE OF THE NOISE (i.e. noise^2)
            if (ct ne 0) then noise2[izero] = 0.0 ;The noise where the weights are zero are set to zero.
            ; and now the same step is done for the two sky nod error frames...
            izero = where(sky1 eq 0.0,ct) ;***
            if (ct ne 0) then sky1[izero] = 1.0 ;***
            skynoise1 = 1.0 / sky1 ;***
            if (ct ne 0) then skynoise1[izero] = 0.0 ;***
            izero = where(sky2 eq 0.0,ct) ;***
            if (ct ne 0) then sky2[izero] = 1.0 ;***
            skynoise2 = 1.0 / sky2 ;***
            if (ct ne 0) then skynoise2[izero] = 0.0 ;***
            izero = 0.0         ;space-saver
            weight = 0.0        ;space-saver
            sky1 = 0.0 ;***
            sky2 = 0.0 ;***

            ptow = 'ptow'+dat[j]
            mask = 'mask'+dat[j]

            wavezp = sxpar(dataheader,'WAVEZP',count=cnt)
            if (cnt eq 0) then wavezp = 0.0 ;if the WAVEZP header element isn't found
            if (wavezp gt 1.0) or (wavezp lt -1.0) then begin
                print,'This wavelength zeropoint is VERY large!'
                print,'Check your fits file!!!'
                print,'File: '+files[j]
                print,wavezp
                stop
            endif
            print,'The wavelength offset is '+strn(wavezp)

            helioc = sxpar(dataheader,'VHELIO',count=cnt)
            if (cnt eq 0) then helioc = 0.0
            print, 'The heliocentric correction is '+strn(helioc)
            
            airmass = sxpar(dataheader,'AIRMASS',count=cnt)
            if (cnt eq 0) then airmass = 1.0 ;if the airmass term isn't written into the header
            print,'The airmass is '+strn(airmass)
        
            readcol,silent=1,dotR[j],format='i,x,x,x,x,x',fiber
            bfi = where(fiber eq -1, count3) ;The bad fibers from PIPE1 are flagged.
            readcol,silent=1,mask,format='i,x',y
            y = y-1             ;the fiber location is now an index

            if (IFU eq'vp2') then tfibers = 246.0 - count3 ;***
            if (IFU eq'vp1') then tfibers = 247.0 - count3 ;***
            tfibers = uint(tfibers) ;***

            if (count3 gt 0) then begin
                for k=0,count3-1 do begin
                    data[*,y[bfi[k]]-step:y[bfi[k]]+step] = -666.0
                    noise2[*,y[bfi[k]]-step:y[bfi[k]]+step] = -666.0 ;these are flagged as -666 flags to go through realignF. they are flagged as = 0.0 later
                    skynoise1[*,y[bfi[k]]-step:y[bfi[k]]+step] = -666.0 ;***
                    skynoise2[*,y[bfi[k]]-step:y[bfi[k]]+step] = -666.0 ;***
                endfor
            endif

            print,'Realigning file '+files[j]
            dataA = realignf(data,ptow,mask,wi,wdisp,aperture,FACTOR=factor,WZP=wavezp,helioC=helioc,VAC=vac) ; the factor is taken into account within realignF. therefore, it wants the wdisp value, NOT wdispF
            n1 = n_elements(dataA[*,0])
            data = 0.0
            datafilenames[j] = files[j]+filesuffix+'A.fits'
            
            if (j eq 0) then begin
                axis1 = n_elements(dataA[*,0])
                axis2 = n_elements(dataA[0,*])
                sxaddpar,headout,'SIMPLE','T'
                sxaddpar,headout,'BITPIX',-64,' Number of bits per data pixel'                  
                sxaddpar,headout,'NAXIS',2,' Number of data axes'
                sxaddpar,headout,'NAXIS1',axis1,' Old size * Factor'
                sxaddpar,headout,'NAXIS2',axis2,' Old size'
                sxaddpar,headout,'EXTEND','T',' FITS data may contain extensions'
                sxaddpar,headout,'CRVAL1',wi,' Initial wavelength'
                sxaddpar,headout,'CDELT1',wdispF,' Wavelength dispersion term'
                if (vac eq 'yes') then sxaddpar,headout,'COMMENT','  The spectra are GIVEN IN VACUUM WAVELENTHS'
                if (vac eq 'no') then sxaddpar,headout,'COMMENT','  The spectra are GIVEN IN AIR WAVELENTHS'
                sxaddpar,headout,'COMMENT','  This spectra has been interpolated to a common wavelength solution'
                sxaddpar,headout,'COMMENT','  This spectra was generated by PIPE2'
                sxaddpar,headout,'COMMENT','  See the original header for other important exposure information'
                datahead = headout
                wave = dindgen(n1) * (wdisp/factor) + wi
            endif

            writefits,datafilenames[j],dataA,datahead
            dataA = 0.0

            noiseA = realignf(noise2,ptow,mask,wi,wdisp,aperture,FACTOR=factor,WZP=wavezp,helioC=helioc)
            noise2 = 0.0
            noisefilenames[j] = files[j]+'perA.fits'

            if (j eq 0) then begin
                sxaddpar,headout,'COMMENT','  This is the interpolated noise^2 frame = 1.0 / pefw'
                weighthead = headout
            endif

            writefits,noisefilenames[j],noiseA,weighthead
            noiseA = 0.0

            skynoise1A = realignf(skynoise1,ptow,mask,wi,wdisp,aperture,FACTOR=factor,WZP=wavezp,helioC=helioc) ;***
            skynoise2A = realignf(skynoise2,ptow,mask,wi,wdisp,aperture,FACTOR=factor,WZP=wavezp,helioC=helioc) ;***
            skyfilenames[0,j] = skyn1[j]+'perA.fits' ;***
            skyfilenames[1,j] = skyn2[j]+'perA.fits' ;***
            writefits,skyfilenames[0,j],skynoise1A,weighthead ;***
            writefits,skyfilenames[1,j],skynoise2A,weighthead ;***
            skynoise1A = 0.0 ;***
            skynoise1A = 0.0 ;***

        endfor

    endif else begin ;ends the linear interpolation part.

        for j=0,n3-1 do begin
            datafilenames[j] = files[j]+filesuffix+'A.fits'
            noisefilenames[j] = files[j]+'perA.fits'
            skyfilenames[0,j] = skyn1[j]+'perA.fits' ;***
            skyfilenames[1,j] = skyn2[j]+'perA.fits' ;***
        endfor

    endelse

;-------------------------------------------------------------------------------
; This is the end of the interpolation section and the beginning of
; the binning section.
;-------------------------------------------------------------------------------

    for j=0,n0-1 do begin      ;a loop though each bin in the bin.list
        binname = binlist[j]
        print,''
        print,'Working on '+binname+'...'
        readcol,silent=1,binname,format='a',fibs
        ikeep = where(fibs ne -1)
        fibs = fibs[ikeep] ;the -1's are tossed from the fiber list, if there are any
        bin = strarr(2,n_elements(fibs))
        for k=0,n_elements(fibs)-1 do bin[*,k] = strsplit(fibs[k],'_',/extract)
        
        cntr = 0
        cntpoint = 0 ;a counter that counts up the number of pointings you have moved through.

        repeat begin
            
            ipp = where(pointing eq pointarr[cntr],countP) ;if the pointing exists
            ibb = where(bin[1,*] eq pointarr[cntr],countB) ;if the pointing is used in a given bin
            
            if (countB gt 0) then begin
                print,strn(countB)+' fiber(s) from pointing '$
                  +pointarr[cntr]+' were found for '+binname
                iff = bin[0,ibb]-1 ;the index of the fibers found for a given bin
            endif
            
            if (countP gt 0) and (countB gt 0) then begin ;if both the pointing exists AND the point is used
                
                for k=0,n_elements(ipp)-1 do begin ;a loop over each pointing (i.e. every 20 min. exp)
                    i1 = ipp[k]
                    readcol,silent=1,dotR[i1],skipline=0,format='i,f,f,f,f,f',fiber,fibrad,dRA,dDEC,lumN,lumR
                    mask = 'mask'+dat[i1]
                    readcol,silent=1,mask,format='i,x',y
                    y = y - 1
                    imm = y[iff] ;the index for the center of ALL THE FIBERS in THIS POINTING
                    
                    imma = intarr(n_elements(imm)*aperture) ;the index for all the ROWS in a pointing
                    cntr2 = 0
                    for l = 0, n_elements(imma)-1, 5 do begin
                        imma[l]   = imm[cntr2]-2
                        imma[l+1] = imm[cntr2]-1
                        imma[l+2] = imm[cntr2]
                        imma[l+3] = imm[cntr2]+1
                        imma[l+4] = imm[cntr2]+2
                        cntr2 = cntr2 + 1
                    endfor
 
                    print,'reading in file '+datafilenames[i1]
                    data = readfits(datafilenames[i1],dataheader,/silent)
                    data = float(data)
                    noise2 = readfits(noisefilenames[i1],/silent)
                    noise2 = float(noise2)
                    sky1 = readfits(skyfilenames[0,i1],/silent) ;***
                    sky2 = readfits(skyfilenames[1,i1],/silent) ;***
                    n1 = n_elements(data[*,0])

                    if (k eq 0) then begin

                        da = data[*,imma]  ;the relevant fibers are selected
                        wa = noise2[*,imma] ;the relevant fibers are selected
                        s1a = sky1[*,imma] ;***
                        s2a = sky2[*,imma] ;***
                        fibnorm = median(da,dim=1) ;the fiber normalization is now done on A SINGLE SCIENCE EXPOSURE.
                        airmass = sxpar(dataheader,'AIRMASS')
                        airmassarr = replicate(airmass,n_elements(imm)*aperture)

                    endif else begin

                        da = [[da],[data[*,imma]]]
                        wa = [[wa],[noise2[*,imma]]]
                        s1a = [[s1a],[sky1[*,imma]]] ;***
                        s2a = [[s2a],[sky2[*,imma]]] ;***
                        fibnorm = [fibnorm,median(data[*,imma],dim=1)]
                        airmass = sxpar(dataheader,'AIRMASS')
                        airmassarr = [airmassarr,replicate(airmass,n_elements(imm)*aperture)]

                    endelse

                    for m=0,n_elements(iff)-1 do begin ;a loop over the index of the fibers found in a given pointing. This is bookkeeping for the dRA and dDec offsets that get written out later.
                        i3 = iff[m]
                        if (fiber[i3] ne -1.0 and m eq 0 and cntpoint eq 0 and k eq 0) then begin
                            alloffsets = [fiber[i3], dRA[i3], dDEC[i3], fibrad[i3], lumR[i3]]
                        endif else begin
                            if (fiber[i3] ne -1.0) then begin
                                piece = [fiber[i3], dRA[i3], dDEC[i3], fibrad[i3], lumR[i3]]
                                if (n_elements(alloffsets) eq 0) then alloffsets = piece else $
                                alloffsets = [[alloffsets],[piece]]
                            endif
                        endelse
                    endfor

                endfor
                
                if (cntpoint eq 0) then begin
                    alldata = da
                    allnoise2 = wa
                    allsky1 = s1a ;***
                    allsky2 = s2a ;***
                    allfibnorm = fibnorm
                    allairmass = airmassarr
                endif else begin
                    alldata = [[alldata],[da]] ;the data and noise get concatenated with the other pointings
                    allnoise2 = [[allnoise2],[wa]]
                    allsky1 = [[allsky1],[s1a]] ;***
                    allsky2 = [[allsky2],[s2a]] ;***
                    allfibnorm = [allfibnorm,fibnorm]
                    allairmass = [allairmass,airmassarr]
                endelse

                cntpoint = cntpoint + 1

            endif
            
            cntr = cntr + 1
            
        endrep until (cntr eq n_elements(pointarr))
        data = 0.0 ;space savers
        noise2 = 0.0
        sky1 = 0.0 ;***
        sky2 = 0.0 ;***

;-------------------------------------------------------------------------------
; The bad fibers (i.e. fibers rejected during pipe1.pro) still exist
; in the arrays. This is a necessary step as the fiber indexing done
; above gets confused if you REMOVE the bad fibers from the data. But,
; now they are just baggage and are tossed from all the arrays.
;-------------------------------------------------------------------------------

        ineg = where(allfibnorm le 0.0,countneg) ;04/22/2011: To catch the possibility of a negative fiber normalization (due to oversubtraction of the sky, most likely)
        if (countneg gt 0) then alldata[*,ineg] = -666 ;04/22/2011: Any oversubtracted fiber row is set to -666 and later rejected

        flags = median(alldata,dimension=1)
        index = where(flags ne -666,count)
        print,''
        print,'The number of total rows in '+binname+' is: '+strn(n_elements(alldata[0,*]))
        print,''
        print,'The number of rows kept, after rejecting bad fibers is: '+strn(count)
        print,''
        print,'The number of rows rejected for over-subtraction is: '+strn(countneg)
        print,''
        if (count eq 0) then begin
            print,'All the data has been rejected for this bin!'
            goto,jumpend
        endif
            
;        if (countneg gt 0) then pause

        nFINAL = count ;the number of kept fiber rows
        alldata = alldata[*,index] ;for large arrays, this step is a killer. 
        allnoise2 = allnoise2[*,index]
        allsky1 = allsky1[*,index] ;***
        allsky2 = allsky2[*,index] ;***
        allfibnorm = allfibnorm[index]
        allairmass = allairmass[index]
        allfibnorm = allfibnorm/max(allfibnorm) ;the weights get normalized to 1.

        if (j eq 0 and skipintpol eq 'yes') then begin
            axis1 = n_elements(alldata[*,0])
            sxaddpar,headout,'SIMPLE','T'
            sxaddpar,headout,'BITPIX',-32,' Number of bits per data pixel'                  
            sxaddpar,headout,'NAXIS',2,' Number of data axes'
            sxaddpar,headout,'NAXIS1',axis1,' Old size * Factor'
            sxaddpar,headout,'NAXIS2',1,' 1-D spectra'
            sxaddpar,headout,'EXTEND','T',' FITS data may contain extensions'
            sxaddpar,headout,'CRVAL1',wi,' Initial wavelength'
            sxaddpar,headout,'CDELT1',wdispF,' Wavelength dispersion term'
            sxaddpar,headout,'COMMENT','  This spectra has been interpolated to a common wavelength solution'
            sxaddpar,headout,'COMMENT','  This spectra was generated by PIPE2'
            sxaddpar,headout,'COMMENT','  See the original file for other important exposure information'
            datahead = headout
            s2nhead = headout
            sxaddpar,datahead,'COMMENT','  The spectra has been normalized to 1, centered around 5000 AA'
            wave = dindgen(n1) * (wdisp/factor) + wi
        endif

        onelisttrim = strsplit(onelist,'.',/extract)
        outnamed = binname+'_'+onelisttrim[0]+'.fits'
        outnamesn = binname+'_'+onelisttrim[0]+'_S2N.fits'

;-----------------------------------------------------------------------------
;                        THE NEW SIGNAL-TO-NOISE CALCULATION
;-----------------------------------------------------------------------------
; The filename_perA.fits files were created earlier. They are the
; interpolated noise^2 frames (i.e. 1/pefw). This new S/N calculation
; now includes shot noise from the two sky nods.

; The noise calculation has the following format:
; N(pixel)^2 = sum_i allnoise2_i + sum_i (t_sci/t_sky)^2 *
; (n_fibers/total_fiber#)^2 allsky1_i + (same sum for allsky2)

        print,'Completing the signal-to-noise calculation...'

        nfibers = n_elements(alldata[0,*])/5.0 ;***
        S2N = fltarr(n1)
        for s=0,n1-1 do begin ;a loop over wavelength
            signal = total(alldata[s,*])
            noise2 = total(allnoise2[s,*]) + ((tsci/tsky)^2 * (nfibers/tfibers)^2 * total(allsky1[s,*])) + $
              ((tsci/tsky)^2 * (nfibers/tfibers)^2 * total(allsky2[s,*]))
            noise = sqrt(noise2)
            S2N[s] = signal / noise
        endfor
        writefits,outnamesn,S2N,s2nhead

;------------------------------------------------------------------------------
;     THE INSTRUMENTAL RESPONSE AND TRANSPARANCY IS BACKED OUT OF THE DATA
;------------------------------------------------------------------------------
        
        if (IR eq 'yes') then begin ; The instrumental response is backed out of the data (added on 12/22/10)
            print,''
            print,'Adjusting '+outnamed+' for instrumental response...'
            if (j eq 0) then begin
                corr1 = interpol(fluxcor,fluxwave,wave)
                corr2 = interpol(transcor,fluxwave,wave)
            endif
            if (TRANS eq 'yes') then begin
                for k=0,nFINAL-1 do alldata[*,k] = alldata[*,k] * corr1 * (10^(0.4 * corr2 *allairmass[k]))
            endif
            if (TRANS eq 'no') then begin
                for k=0,nFINAL-1 do alldata[*,k] = alldata[*,k] * corr1
            endif
        endif
        
; The data is now scaled so as to pass through the biweight in a
; satisfying manner.

        i666 = where(alldata le -666,countbad) ;this index are for cosmic rays. The bad fibers are gone at this point.
        for k=0,nFINAL-1 do alldata[*,k] = alldata[*,k] / allfibnorm[k] ;the data is scaled to be at approximately the same flux.
        if (countbad ne 0) then alldata[i666] = -666.0 ;flags are preserved through the normalization

;------------------------------------------------------------------------------
        if (COMB eq 'biwt') then begin
;                                  THE BIWEIGHT
;------------------------------------------------------------------------------

            if (j eq 0) then sxaddpar,datahead,'COMMENT','  The spectra was made via the biweightF function'
            temp1 = alldata
;            outdata = biweightF(temp1,plotF='on')
;            pause
            outdata = biweightf(temp1,plotF='off')
            temp1 = 0.0
;            alldata = 0.0
        endif

;------------------------------------------------------------------------------
        if (COMB eq 'mean') then begin
;                            THE WEIGHTED AVERAGE
;------------------------------------------------------------------------------

            if (j eq 0) then sxaddpar,datahead,'COMMENT','  The spectra was made via the meanWclipF function'
            temp1 = alldata
            temp2 = 1.0 / allnoise2 ;the meanWclipF routine expects weights, not noise estimates
            outdata = fltarr(n1)
            print,'Running the weighted average...'
            for k=0,n1-1 do outdata[k] = meanWclipF(temp1[k,*],temp2[k,*])
            temp1 = 0.0
            temp2 = 0.0
;            alldata = 0.0
        endif
        
;------------------------------------------------------------------------------
;                       THE FINAL SPECTRA IS NORMALIZED TO 1
;------------------------------------------------------------------------------

        outdata = float(outdata)
        i5000 = where(wave ge 4950.0 and wave le 5050.0)
        norm = median(outdata[i5000])
        outdata = outdata / norm
        writefits,outnamed,outdata,datahead

;------------------------------------------------------------------------------
;                    THE dRA AND dDEC CALCULATION IS MADE
;------------------------------------------------------------------------------

; MODIFIED ON 04/11/09: The old method doesn't work as bins commonly
; combine fibers across the major axis. The new routine now simply
; returns the dRA and dDec values for ALL THE FIBERS in each bin. 

; The routine returns a text file in the following format...
; F #             dRA          dDec          Rad          NFlux

        if (p eq 0) then begin
            free_lun,5
            openw,5,binname+'.coord'
            printf,5,'F #','dRA','dDec','Rad','Med.Flux',$
              format='(a4,2x,11x,a3,10x,a4,10x,a3,7x,a8)'
            for k=0,n_elements(alloffsets[0,*])-1 $
              do printf,5,alloffsets[*,k],$
              format='(i4,2x,f14.5,f14.5,f13.5,f15.5)'
            free_lun,5
        endif

;------------------------------------------------------------------------------
;                THE LAST FIGURES GET PLOTTED, IF TURNED ON.
;------------------------------------------------------------------------------

        if (lastplot eq 'on') then begin
            if (j eq 0) then window,8,retain=2 else wset,8
            plot,wave,S2N,title='The signal-to-noise for '+binname,$
              xrange=[3700,5500],xtitle='Wavelength (A)',/xstyle,$
              ytitle='S/N',/ynozero
            if (j eq 0) then window,10,retain=2 else wset,10
            plot,wave,outdata,title='The final spectra for '+binname,$
              xrange=[3700,5500],xtitle='Wavelength (A)',/xstyle,$
              ytitle='Flux',/ynozero
        endif
        
    jumpend:

    endfor ;index = j (The end of the bin loop.)
endfor ;index = p (The end of the list loop.)

print,''
print,'PIPE2 has finished successfully.'
print,'The next ENTER deletes all open plots.'
print,''

pause
while !d.window ne -1 do wdelete, !d.window

stop
end

;-----------------------------------------------------------------------------
; The old S/N calculation
;-----------------------------------------------------------------------------
        
;        print,'Completing the signal-to-noise calculation...'
;        S2N = fltarr(n1)

;        onelisttrim = strsplit(onelist,'.',/extract)
;        outnamed = binname+'_'+onelisttrim[0]+'.fits'
;        outnamesn = binname+'_'+onelisttrim[0]+'_S2N.fits'
;        outnamesnA = binname+'_'+onelisttrim[0]+'_S2Nalt.fits'

; MODIFICATION MADE ON NOV 30, 2011: There are now 2 different s2n
; values calculated. The first, named S2Nalt, is the classic formula:

;    ***      S/N_tot = sum_i S_i / sqrt (sum_i N_i^2)     ***

; and the second (S2N) is of the form        

;    ***      S/N_tot = sqrt ( sum_i (S_i/N_i)^2)     ***
        
; The second is per Cappellari's 2003 Adaptive binning paper and
; accounts for a weighting of the form w_i = S_i / N_i^2. This is,
; incidentally, the one that has been applied for the past year or so...

; VERSION 1: (S2N.fits)
;        for s=0,n1-1 do begin ;a loop over wavelength
;            temp = dblarr(n_elements(alldata[s,*]))
;            for ss=0,n_elements(alldata[s,*])-1 do begin
;                if (alldata[s,ss] eq -666) then top = 0.0 else top = alldata[s,ss]^2
;                bot = allnoise2[s,ss]
;                temp[ss] = top / bot
;            endfor
;            S2N[s] = sqrt(total(temp))
;        endfor
        
;        S2Nsmooth = smooth(S2N,21,/edge_truncate,/nan)
;        S2N = S2Nsmooth
;        writefits,outnamesn,S2N,s2nhead

; VERSION 2: (S2Nalt.fits)
;        for s=0,n1-1 do begin ;a loop over wavelength
;            temp = dblarr(n_elements(alldata[s,*]))
;            for ss=0,n_elements(alldata[s,*])-1 do begin
;                if (alldata[s,ss] eq -666) then top = 0.0 else top = alldata[s,ss]^2
;                bot = allnoise2[s,ss]
;                temp[ss] = top / bot
;            endfor
;            S2N[s] = sqrt(total(temp))
;        endfor

;        S2Nsmooth = smooth(S2N,21,/edge_truncate,/nan)
;        S2N = S2Nsmooth
;        writefits,outnamesnA,S2N,s2nhead

;*********************************************
;       The very old S/N calculations...        
;*********************************************
        
; This version of the S/N calculation follows M. Cappellari's version
; laid out in a few of his papers on Voronoi binning. It is of the
; form:
;       The signal to noise calculation is now done as follows...
;       S/N_tot = (sum(s_ij/(w_ij*N_ij^2))) / (sqrt(sum(1/(w_ij*N_ij)^2)))

;        for s=0,n1-1 do begin
;            top = total(alldata[s,*] / allnoise2[s,*])
;            t = 1.0 / (allfibnorm^2 * allnoise2[s,*])
;            bottom = sqrt(total(t))
;            S2N[s] = top / bottom
;        endfor
        
;       i = sum over the number of fibers
;       j = sum over the aperture size
;       Yet for me, I am giving EACH ROW it's own weight, so the w's
;       get an index over j.

; NOTE: This is a mess, and should simplify to sqrt ( sum (s/n)^2)

;        onelisttrim = strsplit(onelist,'.',/extract)
;        outnamed = binname+'_'+onelisttrim[0]+'.fits'
;        outnamesn = binname+'_'+onelisttrim[0]+'_s2n.fits'

;        flags = where(alldata eq -666,count)
;        if (count ne 0) then alldata[flags] = 0.0
        
;        signal = dblarr(n1)
;        noise = dblarr(n1)
;        for s=0,n1-1 do signal[s] = total(alldata[s,*])
;        for s=0,n1-1 do begin
;            ind = where(allwgts[s,*] ne 0.0)
;            noise[s] = sqrt(total(1/allwgts[s,ind]))
;        endfor

;        S2N = signal / noise
;        writefits,outnamesn,S2N
;*********************************************        

