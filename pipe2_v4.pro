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
pro pipe2, LIST=listoflists, BIN=yourbinlist,$
              COMB=COMB, IR=IR, TRANS=TRANS,$
              SKIPINTERPOL=skipintpol, VAC=vac
;*****************************************************************

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
wi = 3550.0 ;I have bounced this around in the past, but am now fixing it to 3550, much is about the median of all my data.
wdisp = 1.125 ;(M87:1.125) NOTE: The 'factor' term makes this a 1.125/3.0 = 0.375 dispersion term
factor = 1.0 ;this allows for super-sampling of the wavelength solution by changing the 'wdisp' value to wdisp/factor
aperture = 5 ;The size of the fiber profile, in pixels.
lastplot = 'on' ;This just plots the final spectra for each bin once complete.
filesuffix = 'pefsm' ;This dictates whether you're combining the pefs or pefsm files
pointarr = ['a1','a2','a3','a','b','c','d','e','f','g','h','D1','D2','D3','d1','d2','d3'] ;all the available pointings
swch1 = 'on' & swch2=swch1 & swch3=swch1 & swch4=swch1 & swch5=swch1 & swch6=swch1 & swch7=swch1 & swch8=swch1 & swch9=swch1 & swch10=swch1
COMP = 'areto' ;change to grad78 when running on your grad box (dictates where to find the flux correction files)
tsci = 1200.0 ;***
tsky = 300.0 ;***
normspec = 'no' ;set to 'yes' to normalize the final spectra to have a value of ~1 in the red
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
;The list is composed of 6 columns. The first is the data
;file, sans any pefsm.fits. The second indicates the pointing,
;corresponding to the binlist letters. The third column is the
;corresponding ptow.dat AND mask.dat files. This will be the name you
;gave VACCINE which gets used in the naming convention for the data
;for that run. PIPE2.PRO will append the necessary 'mask' and 'ptow'
;as needed. The 4th column are the "name.R" files that are output from 
;PIPE1.PRO. 

;This is a list of the form:
;  jm0123cb_   a   _jan08_n2.dat    jm0123.R   jm0122_   jm0124_   vp1   old
;  jm0125cd_   a   _jan08_n2.dat    jm0125.R   jm0124_   jm0126_   vp1   old
;  jm3043_     e   _march11_n5.dat  jm3043.R   jm3042_   jm3045_   vp2   new
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

if (n_elements(COMB) eq 0) then COMB = 'biwt' ;the default combination method is the biweight

if (n_elements(listoflists) eq 0) then listoflists = 'data.list'
readcol,silent=1,listoflists,format='A',allists

if (n_elements(skipintpol) eq 0) then skipintpol = 'no'

;The bins are read in
if (n_elements(yourbinlist) eq 0) then readcol,silent=1,'bin.list',format='a',binlist $
else readcol,silent=1,yourbinlist,format='a',binlist
n0 = n_elements(binlist)

for p = 0, n_elements(allists) - 1 do begin ;a loop through each of the lists. THE LINEAR INTERPOLATION FOR WAVELENGTH
    onelist = allists[p]
    print,''
    print,'Working on list '+onelist+'...'
    print,''
    readcol,silent=1,onelist,format='a,a,a,a,a,a,a,a',files,pointing,dat,dotR,skyn1,skyn2,IFUlist,FOCALREDUCERlist

    n3 = n_elements(files) ;the number of files (i.e. number of 20 minute pointings)

    datafilenames = strarr(n3) ;the name of the aligned output files (used for reading in later)
    noisefilenames = strarr(n3) ;same as above
    skyfilenames = strarr(2,n3) ;***

    if (skipintpol eq 'no') then begin ;the file is realigned and written out

        for j=0,n3-1 do begin ;a loop through each file (the length of the datalist)

            print,''
            print,'Reading in file '+files[j]
            data = readfits(files[j]+filesuffix+'.fits',/silent,dataheader)
            weight = readfits(files[j]+'pefw.fits',/silent,noiseheader)
            sky1 = readfits(skyn1[j]+'pefw.fits',/silent,sky1header) ;***
            sky2 = readfits(skyn2[j]+'pefw.fits',/silent,sky2header) ;***

            izero = where(weight eq 0.0,ct)
            if (ct ne 0) then weight[izero] = 1.0
            noise = sqrt(1.0 / weight)
            if (ct ne 0) then noise[izero] = 0.0 ;The noise where the weights are zero are set to zero.
            izero = where(sky1 eq 0.0,ct) ;***
            if (ct ne 0) then sky1[izero] = 1.0 ;***
            skynoise1 = sqrt(1.0 / sky1) ;***
            if (ct ne 0) then skynoise1[izero] = 0.0 ;***
            izero = where(sky2 eq 0.0,ct) ;***
            if (ct ne 0) then sky2[izero] = 1.0 ;***
            skynoise2 = sqrt(1.0 / sky2) ;***
            if (ct ne 0) then skynoise2[izero] = 0.0 ;***
            izero = 0.0         ;space-saver
            weight = 0.0        ;space-saver
            sky1 = 0.0 ;***
            sky2 = 0.0 ;***

            ptow = 'ptow'+dat[j]
            mask = 'mask'+dat[j]

            wavezp = sxpar(dataheader,'WAVEZP',count=cnt)
            if (cnt eq 0) then wavezp = 0.0 ;if the WAVEZP header element isn't found
            if (wavezp gt 2.0) or (wavezp lt -2.0) then begin
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

            if (count3 gt 0) then begin
               for k=0,count3-1 do begin
                  data[*,y[bfi[k]]-step-1:y[bfi[k]]+step+1] = -666.0        ; Bad fibers are flagged as -666 flags to go through realignF. they are flagged as = 0.0 later
                  noise[*,y[bfi[k]]-step-1:y[bfi[k]]+step+1] = -666.0
;                  skynoise1[*,y[bfi[k]]-step-1:y[bfi[k]]+step+1] = -666.0
;                  skynoise2[*,y[bfi[k]]-step-1:y[bfi[k]]+step+1] = -666.0
               endfor
            endif
            
            print,''
            print,'Realigning file '+files[j]
            dataA = realignf(data,ptow,mask,wi,wdisp,aperture,FACTOR=factor,WZP=wavezp,helioC=helioc,VAC=vac) ; the factor is taken into account within realignF. therefore, it wants the wdisp value, NOT wdispF
            n1 = n_elements(dataA[*,0])
            data = 0.0
            
            temp = median(dataA,dim=1)
            ibad = where(temp eq -666.0,count3)
            if (count3 gt 0) then dataA[*,ibad] = -667.0 ;changing the flags to -667 so as to differentiate bad pixels from bad fibers 
            datafilenames[j] = files[j]+filesuffix+'A.fits'

            if (swch1 eq 'on') then begin
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
;                if (IR eq 'yes') then sxaddpar,headout,'COMMENT','  This spectra has been corrected for instrumental response'
;                if (TRANS eq 'yes') then sxaddpar,headout,'COMMENT','  This spectra has been corrected for airmass transparency'
                sxaddpar,headout,'COMMENT','  This spectra has been interpolated to a common wavelength solution'
                sxaddpar,headout,'COMMENT','  This spectra was generated by PIPE2'
                sxaddpar,headout,'COMMENT','  See the original header for other important exposure information'
                datahead = headout
                wave = dindgen(n1) * (wdisp/factor) + wi
                swch1 = 'off'
            endif

            writefits,datafilenames[j],dataA,datahead

            noiseA = realignf(noise,ptow,mask,wi,wdisp,aperture,FACTOR=factor,WZP=wavezp,helioC=helioc)
            noise = 0.0
            noisefilenames[j] = files[j]+'perA.fits'

            if (count3 gt 0) then noiseA[*,ibad] = -667.0 ;the bad FIBERS are flagged in the noise array as well. I do this because I am rejecting the pixels, so I should also be rejecting the noise.

            if (swch2 eq 'on') then begin
                weighthead = headout
                sxaddpar,weighthead,'COMMENT','  This is the interpolated noise frame = sqrt(1.0 / pefw)'
                swch2 = 'off'
            endif

            writefits,noisefilenames[j],noiseA,weighthead

            skynoise1A = realignf(skynoise1,ptow,mask,wi,wdisp,aperture,FACTOR=factor,WZP=wavezp,helioC=helioc) ;***
            skynoise2A = realignf(skynoise2,ptow,mask,wi,wdisp,aperture,FACTOR=factor,WZP=wavezp,helioC=helioc) ;***
            skyfilenames[0,j] = skyn1[j]+'perA.fits' ;***
            skyfilenames[1,j] = skyn2[j]+'perA.fits' ;***
            temp = median(skynoise1A,dim=1)
            ikeep = where(temp ne 0.0)
            col = median(skynoise1A[*,ikeep],dim=2) ;just a single median-collapsed spectra is written out for the sky noise.
            writefits,skyfilenames[0,j],col,weighthead ;***
            temp = median(skynoise2A,dim=1)
            ikeep = where(temp ne 0.0)
            col = median(skynoise2A[*,ikeep],dim=2)
            writefits,skyfilenames[1,j],col,weighthead ;***
            dataA = 0.0 ;space-saving 
            noiseA = 0.0
            skynoise1A = 0.0 ;***
            skynoise2A = 0.0 ;***
            col = 0.0

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

    wgtcntr = 0
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
            
            ipp = where(pointing eq pointarr[cntr],countP) ;if the pointing exists in your data list
            ibb = where(bin[1,*] eq pointarr[cntr],countB) ;if the pointing is used in a given bin
            
            if (countB gt 0) then begin
                print,strn(countB)+' fiber(s) from pointing '$
                  +pointarr[cntr]+' were found for '+binname
                iff = bin[0,ibb]-1 ;the INDEX of the fibers found for a given bin
            endif
            
            if (countP gt 0) and (countB gt 0) then begin ;if both the pointing exists AND the point is used
                for k=0,n_elements(ipp)-1 do begin ;a loop over each exposure for a given pointing
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
                    noise = readfits(noisefilenames[i1],/silent)
                    noise = float(noise)
                    sky1 = readfits(skyfilenames[0,i1],/silent) ;these are now single spectra
                    sky2 = readfits(skyfilenames[1,i1],/silent) 
                    sky1 = float(sky1)
                    sky2 = float(sky2)
                    n1 = n_elements(data[*,0])

                    if (k eq 0) then begin

                        da = data[*,imma]  ;the relevant fibers are selected
                        wa = noise[*,imma] ;the relevant fibers are selected
                        airmass = sxpar(dataheader,'AIRMASS')
                        airmassarr = replicate(airmass,n_elements(imm)*aperture)
                        ifu = IFUlist[i1]
                        focalreducer = FOCALREDUCERlist[i1]
                        ifu = [ifu,ifu,ifu,ifu,ifu]
                        focalreducer = [focalreducer,focalreducer,focalreducer,focalreducer,focalreducer]
                        skyfactor = (tsci/(2.0*tsky) * (n_elements(imma)/1230.0)) ;the sky fudge factor (1230 = 5 * 246)
                        print,''
                        print,'The skyfactor is: '+strn(skyfactor)
                        tempnorm = median(da,dim=1) ;to catch rejected fibers so as not to include them in the skynoise values.
                        itemp = where(tempnorm ne -667.0,tempcount)
                        skynoise1 = sky1 * float(tempcount) * skyfactor ;the multiplication here accounts for the number of rows per fiber and number of fibers coming from a given frame.
                        skynoise2 = sky2 * float(tempcount) * skyfactor
                        skycounter = intarr(n_elements(imma)-1)-1
                        skycounter = [1,skycounter]

                    endif else begin

                        da = [[da],[data[*,imma]]]
                        wa = [[wa],[noise[*,imma]]]
                        airmass = sxpar(dataheader,'AIRMASS')
                        airmassarr = [airmassarr,replicate(airmass,n_elements(imm)*aperture)]
                        ifuT = IFUlist[i1]
                        focalreducerT = FOCALREDUCERlist[i1]
                        ifu = [ifu,ifuT,ifuT,ifuT,ifuT,ifuT]
                        focalreducer = [focalreducer,focalreducerT,focalreducerT,focalreducerT,focalreducerT,focalreducerT]
                        skyfactor = (tsci/(2.0*tsky) * (n_elements(imma)/1230.0))
                        print,''
                        print,'The skyfactor is: '+strn(skyfactor)
                        tempnorm = median(da,dim=1) ;to catch rejected fibers so as not to include them in the skynoise values.
                        itemp = where(tempnorm ne -667.0,tempcount)
                        skynoise1 = [[skynoise1],[sky1 * float(tempcount) * skyfactor]]
                        skynoise2 = [[skynoise2],[sky2 * float(tempcount) * skyfactor]]
                        temp =  intarr(n_elements(imma)-1)-1
                        skycounter = [skycounter,1,temp]

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
                    allnoise = wa
                    allairmass = airmassarr
                    allskynoise1 = skynoise1
                    allskynoise2 = skynoise2
                    allskycounter = skycounter
                endif else begin
                    alldata = [[alldata],[da]] ;the data and noise get concatenated with the other pointings
                    allnoise = [[allnoise],[wa]]
                    allairmass = [allairmass,airmassarr]
                    allskynoise1 = [[allskynoise1],[skynoise1]]
                    allskynoise2 = [[allskynoise2],[skynoise2]]
                    allskycounter = [allskycounter,skycounter]
                endelse

                cntpoint = cntpoint + 1

            endif
            
            cntr = cntr + 1 ;a counter that is used to search through the pointing array and locate the correct pointings that fall into the bin that is being worked on.
            
        endrep until (cntr eq n_elements(pointarr))
        data = 0 ;space savers
        noise = 0
        skynoise1 = allskynoise1
        skynoise2 = allskynoise2
        allskynoise1 = 0
        allskynoise2 = 0
        
;-------------------------------------------------------------------------------
; This ends the selection of fibers that go into a bin.
; The next step is organizational, developing the relevant headers,
; rejecting bad fibers and masked pixels, correctiing for transmission
; and instrumental response, etc.
;-------------------------------------------------------------------------------

        if (swch9 eq 'on' and skipintpol eq 'yes') then begin
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
            if (vac eq 'yes') then sxaddpar,headout,'COMMENT','  The spectra are GIVEN IN VACUUM WAVELENTHS'
            if (vac eq 'no') then sxaddpar,headout,'COMMENT','  The spectra are GIVEN IN AIR WAVELENTHS'
            if (IR eq 'yes') then sxaddpar,headout,'COMMENT','  This spectra has been corrected for instrumental response'
            if (TRANS eq 'yes') then sxaddpar,headout,'COMMENT','  This spectra has been corrected for airmass transparency'
            sxaddpar,headout,'COMMENT','  This spectra was generated by PIPE2'
            sxaddpar,headout,'COMMENT','  See the original file for other important exposure information'
            datahead = headout
            s2nhead = headout
            wave = dindgen(n1) * (wdisp/factor) + wi
            swch9 = 'off'
        endif

        onelisttrim = strsplit(onelist,'.',/extract)
        outnamed = binname+'_'+onelisttrim[0]+'.fits'

        allfibnorm = median(alldata,dim=1) ;array for the scaling of the spectra
        ineg = where(allfibnorm le 0.0,countneg) ;04/22/2011: To catch the possibility of a negative fiber normalization (due to oversubtraction of the sky, most likely).
        if (countneg gt 0) then begin
           allfibnorm[ineg] = -667.0
           alldata[*,ineg] = -667.0 ;04/22/2011: Any oversubtracted fiber row is set to -666 and later rejected
           allnoise[*,ineg] = -667.0
        endif

        GOODindex = where(allfibnorm ne -667.0,countGOOD) ;an index of the GOOD rows
        nFINAL = countGOOD ;the number of kept fiber rows
        nALL = n_elements(alldata[0,*]) ;the number of all rows still left (this includes bad fibers, which get carried through the S/N calculation for accounting reasons.

        print,''
        print,'The number of total rows in '+binname+' is: '+strn(n_elements(alldata[0,*]))
        print,''
        print,'The number of rows kept, after rejecting bad fibers is: '+strn(countGOOD)
        print,''
        print,'The number of rows rejected for over-subtraction is: '+strn(countneg)
        if (countGOOD eq 0) then begin
            print,'All the data has been rejected for this bin!'
            wgtcntr = wgtcntr + 1
            goto,jumpend
        endif

;        test = float(nFINAL/nALL) ;??????????????????????? THIS IS
;        SUPPOSED TO ACCOUNT FOR NOT INCLUDING THE SKY NOISE FROM
;        FIBERS THAT GET REJECTED. 
;        if test lt 1.0 then begin ;if some fibers are tossed the sky noise needs to be adjusted (we don't want to count that noise when we've tossed a fiber)
;           isky = where(allskycounter eq 1,countsky)
;           badfibercorrection = test + (test * 1.0/countsky)
;        endif ;????????????????????????????

;------------------------------------------------------------------------------
;     THE INSTRUMENTAL RESPONSE AND TRANSPARANCY IS BACKED OUT OF THE DATA
;------------------------------------------------------------------------------
 
        i666 = where(alldata eq -666.0,countbad) ;this index are for cosmic rays. The bad fibers are gone at this point.

        if (IR eq 'yes') then begin ; The instrumental response is backed out of the data (added on 12/22/10)
           print,''
           print,'Adjusting '+outnamed+' for instrumental response...'
           cntr3 = 0
           
           for k=0,nALL-1 do begin
              if allfibnorm[k] eq -667.0 then begin
                 print,'Stepping over a bad fiber!'
                 goto, stepoverbadfiber
              endif

              if focalreducer[cntr3] eq 'old' and IFU[cntr3] eq 'vp1' then $
                 readcol,'~/research/flux_calib/fluxfactor_VP1_OFR_FI1_VPH2.txt',$
                         silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
              if focalreducer[cntr3] eq 'old' and IFU[cntr3] eq 'vp2' then $
                 readcol,'~/research/flux_calib/fluxfactor_VP2_OFR_FI1_VPH2.txt',$
                         silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
              if focalreducer[cntr3] eq 'new' and IFU[cntr3] eq 'vp2' then $
                 readcol,'~/research/flux_calib/fluxfactor_VP2_NFR_FI1_VPH2.txt',$
                         silent=1,skipline=4,f='d,d,d',fluxwave,fluxcor,transcor
              corr1 = interpol(fluxcor,fluxwave,wave)

              if (TRANS eq 'no') then begin
                 alldata[*,k] = alldata[*,k] * corr1
                 allnoise[*,k] = allnoise[*,k] * corr1
                 if allskycounter[k] eq 1 then begin
                    skynoise1[*,cntr3] = skynoise1[*,cntr3] * corr1
                    skynoise2[*,cntr3] = skynoise2[*,cntr3] * corr1
                    cntr3 = cntr3 + 1
                 endif
              endif

              if (TRANS eq 'yes') then begin
                 corr2 = interpol(transcor,fluxwave,wave)
                 alldata[*,k] = alldata[*,k] * corr1 * (10^(0.4 * corr2 *allairmass[k]))
                 allnoise[*,k] = allnoise[*,k] * corr1 * (10^(0.4 * corr2 *allairmass[k]))
                 if allskycounter[k] eq 1 then begin
                    skynoise1[*,cntr3] = skynoise1[*,cntr3] * corr1 * (10^(0.4 * corr2 *allairmass[k]))
                    skynoise2[*,cntr3] = skynoise2[*,cntr3] * corr1 * (10^(0.4 * corr2 *allairmass[k]))
                    cntr3 = cntr3 + 1
                 endif
              endif

              stepoverbadfiber:
           endfor
        endif

;-------------------------------------------------------------------------------
; The bad fibers (i.e. fibers rejected during pipe1.pro) still exist
; in the arrays. This is a necessary step as the fiber indexing done
; above gets confused if you REMOVE the bad fibers from the data. But,
; now they are just baggage and are tossed from all the arrays.
;-------------------------------------------------------------------------------

        alldata = alldata[*,GOODindex]
        allnoise = allnoise[*,GOODindex]
        allfibnorm = allfibnorm[GOODindex]
        allairmass = allairmass[GOODindex]
           
;-----------------------------------------------------------------------------
;                        THE SIGNAL-TO-NOISE CALCULATION
;-----------------------------------------------------------------------------
; The filename_perA.fits files were created earlier. They are the
; interpolated noise frames (i.e. sqrt(1/pefw)). This new S/N calculation
; now includes shot noise from the two sky nods.

        print,''
        print,'Completing the signal-to-noise calculation...'
        if (countbad gt 0) then alldata[i666] = 0.0 ;the signal is set to zero so as not to mess with the S2N calculation.
        if (countbad gt 0) then allnoise[i666] = 0.0 ;the noise is set to zero so as not to be included in the noise calculation.
        if (countbad gt 0) then alldataS[i666] = 0.0 ;the signal is set to zero so as not to mess with the S2N calculation.
        if (countbad gt 0) then allnoise[i666] = 0.0 ;the noise is set to zero so as not to be included in the noise calculation.

        S2N1  = fltarr(n1)
        S2N2  = fltarr(n1)

        for s=0,n1-1 do begin          ;a loop over wavelength
           sig1 = total(alldata[s,*])  ;CLASSICAL S/N CALCULATION
           noi1 = sqrt(total(allnoise[s,*]^2)); + total(skynoise1[s,*]^2) + total(skynoise2[s,*]^2))
           noi2 = sqrt(total(allnoise[s,*]^2) + total(skynoise1[s,*]^2) + total(skynoise2[s,*]^2))
           S2N1[s] = sig1/noi1
           S2N2[s] = sig1/noi2
        endfor
        
        for s=0,n1-1 do begin ;a loop over wavelength
           sn = 0.0
           for ss=0,n_elements(alldata[s,*])-1 do begin ;a loop over the number of elements
              signal = alldata[s,ss]
              noise = sqrt(allnoise2[s,ss]) + skyfactor * (sqrt(allsky1[s,ss]) +sqrt(allsky2[s,ss]))
              value = (signal/noise)^2
              sn = [sn,value]
           endfor
           sn = sn[1:*]
           S2Na[s] = sqrt(total(sn))
        endfor

        outarray = fltarr(n1,7) ;the final array: data, S/N1, S/N2 and wavelength
;        outarray = fltarr(n1,3) ;the final array: data, S/N1, S/N2 and wavelength
        outarray[*,4] = S2N1
        outarray[*,5] = S2N2

        if (countbad gt 0) then begin
           alldata[i666] = -666.0 ;now that the S/N is calculated, the masked pixels are set back to -666               
           allnoise[i666] = -666.0
        endif

        alldataS = alldata
        for k=0,nFINAL-1 do alldataS[*,k] = alldata[*,k] / allfibnorm[k] ;the data is scaled to be at approximately the same flux.
        allnoiseS = alldata
        for k=0,nFINAL-1 do alldataS[*,k] = alldata[*,k] / allfibnorm[k] ;the data is scaled to be at approximately the same flux.


;------------------------------------------------------------------------------
        if (COMB eq 'sum') then begin
;                              A STRAIGHT SUMMATION
;------------------------------------------------------------------------------
            if (swch4 eq 'on') then sxaddpar,datahead,'COMMENT','  The spectra was made via a straight summation'
            temp1 = alldata
            outdata = total(temp1,dim=2)
            swch4 = 'off'
        endif
           
;------------------------------------------------------------------------------
        if (COMB eq 'biwt') then begin
;                                  THE BIWEIGHT
;------------------------------------------------------------------------------

           if (swch5 eq 'on') then begin
              sxaddpar,datahead,'COMMENT','  The spectra was made via the biweightf function'
              swch5 = 'off'
           endif
           print,''

           temp1 = alldata
           outdata = biweightf(temp1,plotF='off')
           outarray[*,0] = outdata
           temp1 = alldataS
           outdata = biweightf(temp1,plotF='off')
           outarray[*,2] = outdata

        endif

;------------------------------------------------------------------------------
;        if (COMB eq 'mean') then begin
;                            THE WEIGHTED AVERAGE
;------------------------------------------------------------------------------

;THIS IS CURRENTLY UNDER CONSTRUCTION...

;           if (swch6 eq 'on') then begin
;              sxaddpar,datahead,'COMMENT','  The spectra was made via the meanWclipF function'
;              swch6 = 'off'
;           endif
           temp1 = alldata
           wgts1 = 1.0 / (allnoise^2) ;the classic inverse-variance weighting- this should be equal to the pefw.fits file
;           skywgt1 = 1.0/(sky1^2)
;           skywgt2 = 1.0/(sky2^2)
;           wgts2 = wgts1
;           for s=0,n_elements(alldata[0,*])-1 do begin
;              wgts2[s] = total(alldataCOR[*,s]) / total(1.0/wgts1[*,s] + skyfactor * (1.0/skywgt1 + 1.0/skywgt2))
;           endfor
;           outdata2 = fltarr(n1)
           print,'Running the weighted average...'
           for k=0,n1-1 do outdata[k] = meanwclipf(temp1[k,*],wgts1[k,*])
           outarray[*,1] = outdata
           temp1 = alldataS
           for k=0,n1-1 do outdata[k] = meanwclipf(temp1[k,*],wgts1[k,*])
           outarray[*,3] = outdata
;           array[*,4] = outdata
;           temp1 = alldataCOR
;           for k=0,n1-1 do outdata[k] = meanwclipf(temp1[k,*],wgts2[k,*])
;           array[*,6] = outdata
;           temp1 = 0.0
;           temp2 = 0.0
;        endif
           outarray[*,6] = wave

;------------------------------------------------------------------------------
;                       THE FINAL SPECTRA IS NORMALIZED TO 1
;------------------------------------------------------------------------------

       if (normspec eq 'yes') then begin
          if (swch10 eq 'on') then begin
             sxaddpar,datahead,'COMMENT','  The spectra has been normalized to 1, centered around 5000 AA'
             swch10 = 'off'
          endif

           i5000 = where(wave ge 4950.0 and wave le 5050.0)
           norm = median(outarray[i5000,0],dim=1)
           outarray[*,0] = outarray[*,0] / norm
           norm = median(outarray[i5000,1],dim=1)
           outarray[*,1] = outarray[*,1] / norm
           norm = median(outarray[i5000,2],dim=1)
           outarray[*,2] = outarray[*,2] / norm
           norm = median(outarray[i5000,3],dim=1)
           outarray[*,3] = outarray[*,3] / norm
;           outdata = outdata / norm
        endif

;        outarray[*,0] = outdata
;        outarray[*,3] = wave
       finalarray = fltarr(n1,4) ;scaled biweight, scaled inverse-variance, S/N, wavelength
       finalarray[*,0] = outarray[*,3]
       finalarray[*,1] = outarray[*,4]
       finalarray[*,2] = wave
        writefits,outnamed,finalarray,datahead ;THE FINAL ARRAY IS WRITTEN OUT

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
            if (swch7 eq 'on') then begin
               window,8,retain=2
               swch7 = 'off'
            endif else wset,8
            loadct,0,/silent
;            plot,wave,outarray[*,1],title='The signal-to-noise for '+binname,$
;              xrange=[3700,5500],xtitle='Wavelength (A)',/xstyle,$
;              ytitle='S/N',/ynozero
            plot,wave,outarray[*,0],title='The unscaled data '+binname,$
              xrange=[3500,5800],xtitle='Wavelength (A)',/xstyle,$
              ytitle='S/N',/ynozero
;            wait,0.5
;            loadct,4,/silent
;            oplot,wave,outarray[*,2],color=255
;            pause
            loadct,33
            for jjj=0,n_elements(alldata[0,*])-1 do oplot,wave,alldata[*,jjj],color=jjj
;            pause
            loadct,0,/silent
            oplot,wave,outarray[*,0],thick=2
            oplot,wave,outarray[*,1],color=100,thick=2
            if (swch8 eq 'on') then begin
               window,10,retain=2
               swch8 = 'off'
            endif else wset,10
;            plot,wave,outarray[*,0],title='The final spectra for '+binname,$
;              xrange=[3700,5500],xtitle='Wavelength (A)',/xstyle,$
;              ytitle='Flux',/ynozero
;            loadct,4,/silent
;            oplot,wave,outarray[*,1],color=255
;            loadct,0,/silent
            plot,wave,outarray[*,2],title='The scaled data '+binname,$
              xrange=[3500,5800],xtitle='Wavelength (A)',/xstyle,$
              ytitle='Flux',/ynozero
;            loadct,4,/silent
;            oplot,wave,outarray[*,3],color=255
            loadct,33
            for jjj=0,n_elements(alldatas[0,*])-1 do oplot,wave,alldatas[*,jjj],color=jjj
;            pause
            loadct,0,/silent
            oplot,wave,outarray[*,2],thick=2
            oplot,wave,outarray[*,3],color=100,thick=2
            loadct,0,/silent

;AND NOW THE POSTSCRIPT FILES...

            set_plot,'ps'
            device,file=binname+'unscaled.ps',/color
            loadct,0,/silent
            plot,wave,outarray[*,0],title='The unscaled data '+binname,$
              xrange=[3500,5800],xtitle='Wavelength (A)',/xstyle,$
              ytitle='Flux',/ynozero,xthick=2,ythick=2,charthick=2
            loadct,33
            for jjj=0,n_elements(alldata[0,*])-1 do oplot,wave,alldata[*,jjj],color=jjj
            loadct,0,/silent
            oplot,wave,outarray[*,0],color=50,thick=2
            oplot,wave,outarray[*,1],color=150,thick=2
            device,/close_file
            
            set_plot,'ps'
            device,file=binname+'scaled.ps',/color
            loadct,0
            plot,wave,outarray[*,2],title='The scaled data '+binname,$
              xrange=[3500,5800],xtitle='Wavelength (A)',/xstyle,$
              ytitle='Flux',/ynozero,xthick=2,ythick=2,charthick=2
            loadct,33
            for jjj=0,n_elements(alldatas[0,*])-1 do oplot,wave,alldatas[*,jjj],color=jjj
            loadct,0,/silent
            oplot,wave,outarray[*,2],thick=2
            oplot,wave,outarray[*,3],color=100,thick=2
            device,/close_file

            set_plot,'ps'
            device,file=binname+'S2N.ps',/color
            loadct,0,/silent
            plot,wave,outarray[*,4],title='The signal-2-noise '+binname,$
              xrange=[3500,5800],xtitle='Wavelength (A)',/xstyle,$
              ytitle='S/N',/ynozero,/nodata,xthick=2,ythick=2,charthick=2
            loadct,4
            oplot,wave,outarray[*,4],color=60,thick=2
            oplot,wave,outarray[*,5],color=150,thick=2
            device,/close_file
            set_plot,'x'

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
; Cappellari's calculation, based on Robertson's 1986 paper (S2Nalt1)
;-----------------------------------------------------------------------------

;        skyfactor = (tsci/tsky) * (nfibers/tfibers)
;        skyfactor = 1.0
;        skyfactor = 120.0/246.0
;        for s=0,n1-1 do begin ;a loop over wavelength
;           sn = 0.0
;           for ss=0,n_elements(alldata[s,*])-1 do begin ;a loop over the number of elements
;              signal = alldata[s,ss]
;              noise = sqrt(allnoise2[s,ss]) + skyfactor * (sqrt(allsky1[s,ss]) +sqrt(allsky2[s,ss]))
;              value = (signal/noise)^2
;              sn = [sn,value]
;           endfor
;           sn = sn[1:*]
;           S2Na[s] = sqrt(total(sn))
;        endfor
;        writefits,outnamesnA1,S2Na,s2nhead

;        alldata[i666] = -666 ;the bad pixels are put back...

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
;                if (alldata[s,ss] eq -666) then top = 0.0 else top = alldata[s,ss]^2 ;? is this supposed to be ^2 ?
;                bot = allnoise2[s,ss]
;                temp[ss] = top / bot
;            endfor
;            S2N[s] = sqrt(total(temp))
;        endfor
        
;        S2Nsmooth = smooth(S2N,21,/edge_truncate,/nan)
;        S2N = S2Nsmooth
;        writefits,outnamesnA1,S2N,s2nhead

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
